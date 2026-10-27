# Chargefinetuning

## scripts for: Germline Dependent Charge Gain, Loss, and Reversal during Somatic Hypermutation in Antibodies

Vera A. Spanke¹, Katharina B. Kroell¹, Adrian I. Senn¹, Alexander Bujotzek², Guy Georges², Klaus R. Liedl¹

¹: Department of General, Inorganic and Theoretical Chemistry, University of Innsbruck, Innsbruck, Austria
²: Roche Pharma Research and Early Development, Large Molecule Research, Roche Innovation Center Munich, Penzberg, Germany 

As soon as the manuscript is accepted the abstract will be here. Creates Figure 5 of manuscript

## Requirements

- Python 3.9.2
- Packages listed in `requirements.txt`

Install the Python dependencies:
```bash
python -m pip install -r requirements.txt
```

how to use the immcantation docker: 
```bash
docker pull immcantation/suite:4.3.0
docker run -it --workdir /data -v $(pwd):/data:z immcantation/suite:4.3.0 bash 
```

## Additional Input data

prior running the scripts, the .airr version of the paired AbNGS must be downloaded (https://naturalantibody.com/paired-ab-ngs/) and saved in a separated directory

## Usage

Run the relevant script from the project directory in the following order:

1. prepare_all_chains_to_fasta.ipynb
2. in docker: 
```bash
bash heavy_chain_clonotypes.sh 
```
3. repeat for lambda and kappa light chains, change PREFIX in script
4. separate_by_germline_study.ipynb
5. in docker: 
```bash
bash all_observed_expected.sh
```
6. enriched_depleted.ipynb

## Output

The scripts generate analysis results and figure 5 of the manuscript: IMGT position, germline V gene, and study wise observed and expected charge related mutation frequencies


