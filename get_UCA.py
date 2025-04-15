import argparse
import subprocess
import sys
from Bio.Seq import Seq
from multiprocessing import Pool
import numpy as np
import olga.load_model as load_model
import olga.generation_probability as pgen
import pandas as pd 
import random
import os

def read_starting_and_ending_points(file_path):
    with open(file_path, 'r') as file:
        content = file.read().strip()
        starting_point, ending_point = map(int, content.split())
    return starting_point, ending_point

def generate_codon_list(codons, starting_germline, starting_point, end_point, cdr3_only):
    if cdr3_only:
        codons_cdr3 = [(i, starting_germline[i:i+3]) for i in range(starting_point, end_point, 3)]
        positions = [pos for pos, _ in codons_cdr3]
        new_order = random.sample(positions, len(positions))
    else:
        new_order = random.sample(range(len(codons)), len(codons))
    return new_order

def propose_new_combination(current_germline, codon_number, igphyml_df):
    current_codon = current_germline[codon_number:codon_number+3]
    sub_df = igphyml_df[igphyml_df["site"] == codon_number/3]
    sub_df = sub_df[sub_df["codon"] != current_codon]
    new_codon = sub_df.sample(n=1)
    new_germline = current_germline[:codon_number] + new_codon["codon"].values[0] + current_germline[codon_number+3:]
    return new_germline 
    
def get_new_lhood(proposed_combination, starting_point, ending_point, pgen_model, igphyml_df):
    cdr3 = proposed_combination[starting_point:ending_point]
    new_pgen = pgen_model.compute_nt_CDR3_pgen(cdr3)
    if new_pgen == 0:
        new_pgen = new_pgen + 1e-323
    new_pgen = np.log(new_pgen)
    new_codons = [proposed_combination[i:i+3] for i in range(0, len(proposed_combination), 3)]
    parital_likelihoods = []
    for i in range(0, igphyml_df['site'].max()):
        matched_value = igphyml_df.loc[(igphyml_df['site'] == i) & (igphyml_df['codon'] == new_codons[i]), 'partial_likelihood'].values[0]
        parital_likelihoods.append(matched_value)
    new_tree_likelihood = sum(parital_likelihoods)
    new_lhood = new_pgen + new_tree_likelihood
    lhood_vector = [new_pgen, new_tree_likelihood, new_lhood]
    return lhood_vector

def translate_to_amino_acid(codon):
    return str(Seq(codon).translate())

def process_option(option, filtered_df, codon, current_germline, starting_point, ending_point, igphyml_df_new, pgen_model):
    testing_df = filtered_df.iloc[option]
    testing_germline = current_germline 
    testing_germline = [testing_germline[i:i+3] for i in range(0, len(testing_germline), 3)]
    testing_germline[codon // 3] = testing_df['codon']
    testing_germline = ''.join(testing_germline)
    new_lhoods = get_new_lhood(testing_germline, starting_point, ending_point, pgen_model, igphyml_df_new)
    return testing_germline, new_lhoods

def process_option_wrapper(args):
    return process_option(*args)

def get_updated_germline(starting_germline, starting_cdr3, igphyml_df, pgen_model, starting_point, ending_point, cdr3_only=True, max_iter=100, nproc = 1):
    current_codons = [starting_germline[i:i+3] for i in range(0, len(starting_germline), 3)]
    current_germline = starting_germline
    current_lhoods = get_new_lhood(current_germline, starting_point, ending_point, pgen_model, igphyml_df)
    starting_site = starting_point/3
    ending_site = ending_point/3 - 1
    starting_options = igphyml_df[igphyml_df["site"] == starting_site]
    starting_options = starting_options[starting_options['codon'].apply(translate_to_amino_acid) == 'C']
    ending_options = igphyml_df[igphyml_df["site"] == ending_site]
    ending_options = ending_options[ending_options['codon'].apply(translate_to_amino_acid).isin(['W', 'F'])]
    filtered_igphyml_df = igphyml_df[~igphyml_df['site'].isin([starting_site, ending_site])]
    igphyml_df_new = pd.concat([filtered_igphyml_df, starting_options, ending_options])
    igphyml_df_new = igphyml_df_new.sort_values(by='site', ascending=True)
    still_improving = True
    tested_combinations = []
    tested_lhoods_pgen = []
    tested_lhoods_tree = []
    tested_lhoods_combo = []
    tested_iterations = []
    iteration = 0

    while still_improving and iteration < max_iter:
        still_improving = False
        codon_list = generate_codon_list(current_codons, current_germline, starting_point, ending_point, cdr3_only)
        iteration += 1
        for codon in codon_list:           
            filtered_df = igphyml_df_new[igphyml_df_new["site"] == codon/3]
            with Pool(nproc) as pool:
                results = pool.map(process_option_wrapper, [(option, filtered_df, codon, current_germline, starting_point, ending_point, igphyml_df_new, pgen_model) for option in range(0, len(filtered_df))])
            results_array = np.array([result[1] for result in results])
            best_index = np.argmax(results_array[:, 2])
            tested_combinations.append(results[best_index][0])
            tested_lhoods_pgen.append(results[best_index][1][0])
            tested_lhoods_tree.append(results[best_index][1][1])
            tested_lhoods_combo.append(results[best_index][1][2])
            indices = [i for i, x in enumerate(codon_list) if x == codon]
            tested_iterations.extend([f"{iteration}_{index}" for index in indices])
            if results_array[best_index][2] > current_lhoods[2]:
                current_lhoods = list(results_array[best_index])
                current_germline = results[best_index][0]
                still_improving = True

    data = {
    'germline': tested_combinations,
    'pgen_lhood': tested_lhoods_pgen,
    'tree_lhood': tested_lhoods_tree,
    'combo_lhood': tested_lhoods_combo,
    'iteration': tested_iterations,
    }
    data = pd.DataFrame(data)
    return current_germline, current_lhoods, data

if __name__ == '__main__':
    # Initialize the argument parser
    parser = argparse.ArgumentParser(description='Process UCA arguments.')
    parser.add_argument('--clone_ids', required=True, help='A comma-separated string of the clone IDs to get UCA for')
    parser.add_argument('--directory', required=True, help='Directory where the clone data is stored')
    parser.add_argument('--max_iters', type=int, required=True, help='Max number of iterations to run')
    parser.add_argument('--nproc', type=int, required=True, help='Number of processors to use')
    parser.add_argument('--id', required=True, help='The name of the folder that is created to store the data')
    parser.add_argument('--model_folder', required=True, help='The file path to the OLGA model files for IGH')
    parser.add_argument('--model_folder_igk', required=True, help='The file path to the OLGA model files for IGK')
    parser.add_argument('--model_folder_igl', required=True, help='The file path to the OLGA model files for IGL')
    parser.add_argument('--quiet', type=int, default=0, help='Whether to print out the progress/messages of the script')
    parser.add_argument("--starting_germlines", required=True, help="A comma-separated string of the file paths to the starting germline. "
                                                                    "If you have run getTreesAndUCA this file will be named 'olga_testing_germline.txt'")
    parser.add_argument("--junction_locations", required=True, help="A comma-separated string of the file paths to the file containing the starting and ending site number of the junction. "
                                                                    "If you have run getTreesAndUCA this file will be named 'olga_junction_positions.txt'")
    parser.add_argument("--tree_tables", required=True, help="A comma-separated string of the file paths to the tree tables associated with these clones. "
                                                            "If you have run getTreesAndUCA this file will end in '_pars_hlp_rootprobs.txt' or if running with both heavy and light chains it will be 'heavy_table.txt' "
                                                            "'light_table.txt'")
    parser.add_argument("--chains", required=True, help="A comma-separated string of the locus chains (IGH, IGK, or IGL) associated with these clones.")
    args = parser.parse_args()

    if args.quiet > 0:
        print("Arguments: ", args)

    clone_ids_list = args.clone_ids.split(',')
    starting_germlines_list = args.starting_germlines.split(',')
    junction_locations_list = args.junction_locations.split(',')
    tree_table_list = args.tree_tables.split(',')
    chains_list = args.chains.split(',')

    index_table = pd.DataFrame({
        'clone_ids': clone_ids_list,
        'starting_germline': starting_germlines_list,
        'junction_locations': junction_locations_list,
        'tree_table': tree_table_list,
        'chains': chains_list
    })

    for index, row in index_table.iterrows():
        clone_number = row['clone_ids']
        base_string = args.directory + "/" + args.id + "_" + clone_number
        if args.quiet > 0:
            print("\nRunning on clone", clone_number, "with saved data in", base_string)
        
        if not os.path.exists(base_string):
            os.makedirs(base_string)
        chain = row['chains']

        table_path = row['tree_table']    
        igphyml_df = pd.read_csv(table_path, sep = "\t", header = None)
        igphyml_df.columns = ["site", "codon", "partial_likelihood", "nope", "nada", "no", "equilibrium"]
        igphyml_df['value'] = igphyml_df['partial_likelihood'] + np.log(igphyml_df['equilibrium'])

        germline_string = row['starting_germline']
        junction_string = row['junction_locations']
        starting_point, ending_point = read_starting_and_ending_points(junction_string)
        with open(germline_string, "r") as f:
            starting_germline = f.read().strip()
        starting_cdr3 = starting_germline[starting_point:ending_point]

        v_gene = starting_germline[0:starting_point]
        j_gene = starting_germline[ending_point:len(starting_germline)]
        if "N" in v_gene or "N" in j_gene:
            if(args.quiet > 0):
                print("N in the J gene -- adding to igphyml df")
            codons = [starting_germline[i:i+3] for i in range(0, len(starting_germline), 3)]
            indices_with_N = [index for index, value in enumerate(codons) if 'N' in value and not (len(v_gene) // 3 + 1 <= index <= len(codons) - len(j_gene) // 3  - 1)]
            new_rows = []
            for index in indices_with_N:
                sub = igphyml_df[igphyml_df['site'] == index]
                sub_sum = sub['value'].sum()
                new_rows.append({'site': index, 'codon': codons[index], 'partial_likelihood': sub_sum, 'nope': 0, 'nada': 0, 'no': 0, 'equilibrium': 0, 'value': 0})
            igphyml_df = pd.concat([igphyml_df, pd.DataFrame(new_rows)], ignore_index=True)
            igphyml_df = igphyml_df.sort_values(by='site', ascending=True)
            # replace the Ns in the cdr3 with C and stitch it back together
            cdr3 = starting_germline.split(v_gene)
            cdr3  = cdr3[1].split(j_gene)[0]
            starting_germline = v_gene + cdr3.replace("N", "C") + j_gene

        if "N" in starting_cdr3:
            if(args.quiet > 0):
                print("N in starting CDR3 -- replacing with C")
            starting_cdr3 = starting_cdr3.replace("N", "C")
            # replace the Ns in the cdr3 with C and stitch it back together
            starting_germline = v_gene + starting_cdr3 + j_gene

        # split the germline into a list of codons
        codons = [starting_germline[i:i+3] for i in range(0, len(starting_germline), 3)]

        if chain == "IGH":
            model_folder = args.model_folder
        elif chain == "IGK":
            model_folder = args.model_folder_igk
        elif chain == "IGL":
            model_folder = args.model_folder_igl
        else:
            print("Invalid chain type. Please use IGH, IGK, or IGL.")
            sys.exit(1)

        params_file_name = model_folder + '/model_params.txt'
        marginals_file_name = model_folder + '/model_marginals.txt'
        V_anchor_pos_file = model_folder + '/V_gene_CDR3_anchors.csv'
        J_anchor_pos_file = model_folder + '/J_gene_CDR3_anchors.csv'

        if chain == "IGH":
            genomic_data = load_model.GenomicDataVDJ()
            genomic_data.load_igor_genomic_data(params_file_name, V_anchor_pos_file, J_anchor_pos_file)
            generative_model = load_model.GenerativeModelVDJ()
            generative_model.load_and_process_igor_model(marginals_file_name)
            pgen_model = pgen.GenerationProbabilityVDJ(generative_model, genomic_data)
        else:
            genomic_data = load_model.GenomicDataVJ()
            genomic_data.load_igor_genomic_data(params_file_name, V_anchor_pos_file, J_anchor_pos_file)
            generative_model = load_model.GenerativeModelVJ()
            generative_model.load_and_process_igor_model(marginals_file_name)
            pgen_model = pgen.GenerationProbabilityVJ(generative_model, genomic_data)
        

        # TODO: if the model (pgen) files aren't included and need to be read in, the argv calls here will need to be updated
        values = get_updated_germline(starting_germline, starting_cdr3, igphyml_df, pgen_model, starting_point, ending_point, cdr3_only=True, max_iter=int(args.max_iters), nproc = int(args.nproc))

        new_germline = values[0]
        new_lhoods = values[1]
        data = values[2]

        try:
            if chain == "IGH":
                with open(base_string + "/UCA.txt", "w") as f:
                    f.write(new_germline + "\n")
                with open(base_string + "/UCA_lhoods.txt", "w") as f:
                    f.write(str(new_lhoods) + "\n")
                data.to_csv(base_string + "/UCA_data.csv", index=False)
            else:
                with open(base_string + "/UCA_light.txt", "w") as f:
                    f.write(new_germline + "\n")
                with open(base_string + "/UCA_lhoods_light.txt", "w") as f:
                    f.write(str(new_lhoods) + "\n")
                data.to_csv(base_string + "/UCA_data_light.csv", index=False)
        except Exception as e:
            print(f"Error writing output files: {e}")