#/bin/bash
set -ex;

python Simulation_TwoCelltypes.py                     --out-root simulations_2celltypes                     --omega-shrink 0.9 --max-hsq 0.8                     --snps 1000 --metric type_1_error --parallel --method blue --trun-corr --ieqtl 

python Simulation_TwoCelltypes.py                     --out-root simulations_2celltypes                     --omega-shrink 0.9 --max-hsq 0.8                     --snps 1000 --metric type_1_error --parallel --method tram 

python Simulation_TwoCelltypes.py                     --out-root simulations_2celltypes                     --omega-shrink 0.9 --max-hsq 0.8                     --snps 1000 --metric power --parallel --method blue --trun-corr --ieqtl 

python Simulation_TwoCelltypes.py                     --out-root simulations_2celltypes                     --omega-shrink 0.9 --max-hsq 0.8                     --snps 1000 --metric power --parallel --method tram 


python Simulation_MultiCelltypes.py                     --out-root simulations_6celltypes                     --omega-shrink 0.9 --max-hsq 0.8                     --snps 1000 --metric type_1_error --parallel --method blue --trun-corr --ieqtl 

python Simulation_MultiCelltypes.py                     --out-root simulations_6celltypes                     --omega-shrink 0.9 --max-hsq 0.8                     --snps 1000 --metric type_1_error --parallel --method tram


python Simulation_MultiCelltypes.py                     --out-root simulations_6celltypes                     --omega-shrink 0.9 --max-hsq 0.8                     --snps 1000 --metric power --parallel --method blue --trun-corr --ieqtl 

python Simulation_MultiCelltypes.py                     --out-root simulations_6celltypes                     --omega-shrink 0.9 --max-hsq 0.8                     --snps 1000 --metric power --parallel --method tram