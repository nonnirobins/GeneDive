#!/bin/bash

# Initialize variables
ASSEMBLY_DIR=""
HMM_DIR=""
OUTPUT_DIR=""

# Function to show usage
usage() {
    echo "Usage: $0 -a <assembly_directory> -h <hmm_directory> -o <output_directory>"
    exit 1
}

# Parse command-line options
while getopts ":a:h:o:" opt; do
  case $opt in
    a)
      ASSEMBLY_DIR=$OPTARG
      ;;
    h)
      HMM_DIR=$OPTARG
      ;;
    o)
      OUTPUT_DIR=$OPTARG
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      usage
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      usage
      ;;
  esac
done

# Check if all required inputs are provided
if [ -z "$ASSEMBLY_DIR" ] || [ -z "$HMM_DIR" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "All options -a, -h, and -o are required."
    usage
fi

# Function to check if Conda is installed
check_conda_installed () {
    if ! type "conda" > /dev/null; then
        echo "Conda is not installed. Installing Conda..."
        wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O Miniconda.sh
        bash Miniconda.sh -b
        eval "$(/home/$USER/miniconda3/bin/conda shell.bash hook)"
    else
        echo "Conda is already installed."
    fi
}

# Function to create a temporary Conda environment and install dependencies
create_conda_environment_from_yaml () {
    local env_name="genedive_env"
    local yml_file="geneDive_env.yml"

    if ! conda info --envs | grep -q "$env_name"; then
        echo "Creating Conda environment $env_name from $yml_file..."
        conda env create -f "$yml_file"
    else
        echo "Conda environment $env_name already exists."
    fi
    conda activate "$env_name"
}

# Function to remove the temporary Conda environment
cleanup_environment () {
    conda deactivate
}

# Function to run a command and exit on failure
run_command() {
    "$@"
    if [ $? -ne 0 ]; then
        echo "Error executing $1"
        exit 1
    fi
}

# Find Conda base directory and initialize it
CONDA_BASE=$(conda info --base)
if [ -f "$CONDA_BASE/etc/profile.d/conda.sh" ]; then
    source "$CONDA_BASE/etc/profile.d/conda.sh"
else
    echo "Error: unable to initialize Conda."
    exit 1
fi

# Check if Conda is installed and install if necessary
check_conda_installed

# Create a temporary environment and install dependencies
create_conda_environment_from_yaml

# Main processing logic
declare -A pd_table  # Associative array for storing PD values

for assembly_file in "$ASSEMBLY_DIR"/*; do
  for hmm_file in "$HMM_DIR"/*; do
    # Set variables for this iteration
    ASSEMBLY=$assembly_file
    HMM_FILE=$hmm_file
    OUTPUT=$OUTPUT_DIR/$(basename "$assembly_file")_$(basename "$hmm_file")

    # Ensure output directory exists
    mkdir -p $OUTPUT

    # GeneDive analysis
    echo "Running analysis for $(basename "$assembly_file") and $(basename "$hmm_file")"

    # Genome Annotation with Prokka
    echo "Running FragGeneScan for ORF identification..."
    run_command FragGeneScan -s $ASSEMBLY -o $OUTPUT.annotation -w 1 -t complete 

    # HMM Search
    echo "Running HMM search..."
    run_command hmmsearch -A $OUTPUT.hmmsearch_output.sto "$HMM_FILE" $OUTPUT.annotation.faa

    esl-alimask --rf-is-mask $OUTPUT.hmmsearch_output.sto |\
        esl-reformat --mingap afa -|\
        esl-alimanip --lmin 100 - |\
        esl-reformat afa -|\
        cut -d ' ' -f 1 >  $OUTPUT.hmmsearch_output.afa
    cat $OUTPUT.hmmsearch_output.afa

    # Construct Phylogenetic Tree with FastTree
    echo "Constructing phylogenetic tree with FastTree..."
    run_command fasttree $OUTPUT.hmmsearch_output.afa > $OUTPUT.tree_output

    # Python script embedded in bash script
    echo "Calculating Faith's PD..."

    read -r -d '' PYTHON_SCRIPT << 'EOF'
    #!/usr/bin/env python

    import sys
    from ete3 import Tree

    def install_ete3():
        try:
            import ete3
        except ImportError:
            print("Installing ete3...")
            import subprocess
            subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'ete3'])

    def calculate_faith_pd(tree_file):
        # Load the tree
        tree = Tree(tree_file)

        # Calculate Faith's PD based on branch lengths
        faith_pd = sum(node.dist for node in tree.traverse())

        return faith_pd

    def main():
        # Install ete3 if not already installed
        install_ete3()

        # Parse arguments
        if len(sys.argv) != 2:
            print("Usage: python script.py <tree_file>")
            sys.exit(1)

        tree_file = sys.argv[1]

        # Calculate and print Faith's PD
        faith_pd = calculate_faith_pd(tree_file)
        print(faith_pd)

    if __name__ == "__main__":
        main()
EOF


    # Calculate Faith's PD
    faith_pd=$(python -c "$PYTHON_SCRIPT" "$OUTPUT.tree_output")
    pd_table[$(basename "$assembly_file")][$(basename "$hmm_file")]=$faith_pd
  done
done


# Output table generation
echo "Sample,Gene1,Gene2,...,GeneN" > "$OUTPUT_DIR/output_table.csv"
for sample in "${!pd_table[@]}"; do
  echo -n "$sample" >> "$OUTPUT_DIR/output_table.csv"
  for gene in "${!pd_table[$sample][@]}"; do
    echo -n ",${pd_table[$sample][$gene]}" >> "$OUTPUT_DIR/output_table.csv"
  done
  echo >> "$OUTPUT_DIR/output_table.csv"
done

echo "GeneDive analysis complete"
