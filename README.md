# Novel Model Slection Criterion for Inference of Ising Models
## Directory:
1. fi0201.m
	- Primary script for the full forwad Ising and inference problem
	- Calls to supporting subscripts in Mike_Code_4 subfolder  
	- Responsible for defining regime parameters,
	- Intakes the job name and assigned regimes

2. Mike_Code_4
    - Contains the supporting scripts and matlab functions
    - TCS4.m
    	- Forward Ising process, generates the Adjacency matrix and Ising parameters
    - Met_Hast_norm.m
    	- Metropolis Hastings Algo, generates the simulated observation samples
    - PBLLH6.m
    	- Inverse solution, implements the pseudo-likelihood and selection criteria
    	- Also calculates the performance of the methods

3. SWWC0201.slurm
	- SLURM scheduler job submission
	- Assigns the IDs of the regime parameters to be run
