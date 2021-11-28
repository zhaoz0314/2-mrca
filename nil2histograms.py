import msprime
from pathlib import Path

# set parameters
trial_number= 90000 # num of trials
v_bin_number= 200 # num of bins
sample_size=2 # num of individuals
v_ploidy=1 # num of chromosomes of each type
combined_parameters= [0, 1/2, 1, 2, 4, 8, 16, 32, 40, 48, 56, 64] # 2*sequence_length*recombination_rate*population_size
sequence_lengths= [100*10**idx for idx in range(3)] # sequence length
population_sizes= [1000*10**idx for idx in range(3)] # population size
parameter_packages= [[combined_parameter, v_sequence_length, v_population_size,
                       combined_parameter/2/v_sequence_length/v_population_size]
                      for combined_parameter in combined_parameters
                      for v_sequence_length in sequence_lengths
                      for v_population_size in population_sizes] # num of recombinations per base per gen
# relevant_recombination_rate=2*sequence_length*recombination_rate

# define the function that finds the least coalescent time
def mrca_time_api(tree_sequence):
    ca_times= []
    for tree in tree_sequence.trees():
        ca_times.append(tree_sequence.tables.nodes[tree.root].time)
    return(min(ca_times))

# define the function that transforms lists to histograms
def list_2_histogram(f_list,bin_number):
    bin_size= max(f_list)/bin_number
    bin_locations= [(idx+0.5)*bin_size for idx in range(bin_number)]
    bin_heights= [0 for idx in range(bin_number)]
    for entry in f_list:
        bin_idx= int(entry/bin_size)
        if bin_idx== bin_number:
            bin_idx= bin_idx-1
        bin_heights[bin_idx]= bin_heights[bin_idx]+1/len(f_list)/bin_size
    return [bin_locations, bin_heights]

# save the parameter packages
subdirectory= ("{:.1e}".format(trial_number)+ "_"+
           "{:.1e}".format(v_bin_number)+ "/")
Path(subdirectory).mkdir()
with open(subdirectory+ "parameter_packages.txt", 'w') as txtfile:
    txtfile.write(str(parameter_packages))
with open(subdirectory+ "combined_parameters.txt", 'w') as txtfile:
    txtfile.write(str(combined_parameters))
with open(subdirectory+ "sequence_lengths.txt", 'w') as txtfile:
    txtfile.write(str(sequence_lengths))
with open(subdirectory+ "population_sizes.txt", 'w') as txtfile:
    txtfile.write(str(population_sizes))

# run the simulation trial_number times for each parameter package
# and save each histogram
for parameter_package in parameter_packages:
    mrca_times= []
    for idx in range(trial_number):
        tree_sequence_simulated= msprime.sim_ancestry(
            samples= sample_size,
            ploidy= v_ploidy,
            recombination_rate= parameter_package[3],
            sequence_length= parameter_package[1],
            population_size= parameter_package[2],
            record_full_arg= True)
        mrca_times.append(mrca_time_api(tree_sequence_simulated)/parameter_package[2])
    histogram= list_2_histogram(mrca_times, v_bin_number)
    # with open("parameter_package"+datetime.now().strftime("%H_%M_%S")+".txt", 'w') as txtfile:
    filename= ("{:.1e}".format(parameter_package[0])+ "_"+
               "{:.1e}".format(parameter_package[1])+ "_"+
               "{:.1e}".format(parameter_package[2])+ ".txt")
    with open(subdirectory+ filename, 'w') as txtfile:
        txtfile.write(str(histogram))
