import sys
import json
import argparse

def import_json_file(file_path):
    """
    Import JSON data from a file.

    Args:
        file_path (str): Path to the input JSON file.

    Returns:
        dict: JSON data as a dictionary.
    """
    with open(file_path, 'r') as file:
        data = json.load(file)
    return data

def export_json_file(data, file_path):
    """
    Export JSON data to a file.

    Args:
        data (dict): JSON data to be exported.
        file_path (str): Path to the output JSON file.
    """
    with open(file_path, 'w') as file:
        json.dump(data, file, indent=4)

def main():
    """
    Change the mmseqs --min-seq-id parameter in the EcoPhylo Snakemake JSON config file.
    """
    parser = argparse.ArgumentParser(description='Change the mmseqs --min-seq-id parameter in the EcoPhylo Snakemake JSON config file.')

    # Add command-line arguments
    parser.add_argument('--input-config-json', type=str, help='Path to the input JSON file.')
    parser.add_argument('--output-config-json', type=str, help='Path to the output JSON file.')
    parser.add_argument('--min-seq-id', type=int, help='New value for --min-seq-id.')
    parser.add_argument('--cov-mode', type=int, help='New value for mmseq --cov-mode.')
    parser.add_argument('--hmm-list', type=str, help='hmm_list.txt file path')
    parser.add_argument('--SCG', type=str, help='SCG name')
    parser.add_argument('--cami-dataset', type=str, help='CAMI dataset name')
    parser.add_argument('--external-genomes', type=bool, help='Will there be an external-genomes.txt file?')

    # Parse the command-line arguments
    args = parser.parse_args()

    # Import JSON file
    file_path = args.input_config_json
    hmm_list = args.hmm_list
    file = import_json_file(file_path)

    # Modify the parameters in the JSON file
    if args.min_seq_id == 1:
        min_seq_id = 1
    else:
        min_seq_id = args.min_seq_id / 100
    file['cluster_X_percent_sim_mmseqs']['--min-seq-id'] = min_seq_id
    file['cluster_X_percent_sim_mmseqs']['--cov-mode'] = args.cov_mode
    file['cluster_X_percent_sim_mmseqs']['threads'] = 20
    file['hmm_list'] = hmm_list
    file['run_genomes_sanity_check'] = False
   
    samples_txt = "samples_" + args.cami_dataset + ".txt"
    file['samples_txt'] = samples_txt
    metagenomes_txt = "metagenomes_" + args.cami_dataset + ".txt"
    file['metagenomes'] = metagenomes_txt

    if args.external_genomes:
        external_genomes_txt = "external_genomes_" + args.cami_dataset + ".txt"
        file['external_genomes'] = external_genomes_txt
    else:
        file['external_genomes'] = ""

    # Modify EcoPhylo directory structure names

    for key in file['output_dirs']:
        if 'ECOPHYLO_WORKFLOW' in file['output_dirs'][key]:
            if not args.external_genomes:
                file['output_dirs'][key] = file['output_dirs'][key].replace('ECOPHYLO_WORKFLOW', f'ECOPHYLO_WORKFLOW_{args.min_seq_id}_PERCENT_SIM_{args.SCG}_{args.cami_dataset}_NO_EXTERNAL_GENOMES')
            else:
                file['output_dirs'][key] = file['output_dirs'][key].replace('ECOPHYLO_WORKFLOW', f'ECOPHYLO_WORKFLOW_{args.min_seq_id}_PERCENT_SIM_{args.SCG}_{args.cami_dataset}')

    # Modify metagenomics workflow parameters
    file['run_metagenomics_workflow']['clusterize'] = True
    file['run_metagenomics_workflow']['clusterize_submission_params'] = '--exclude \'\''
    file['run_metagenomics_workflow']['snakemake_additional_params'] = '--jobs=80 --resources nodes=80'

    # Export modified JSON data to a file
    output_file_path = args.output_config_json
    export_json_file(file, output_file_path)

if __name__ == '__main__':
    main()
