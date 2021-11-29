import msprime
import matplotlib.pyplot as plt

trial_number= 100000
v_bin_number= 200
sample_size= 2
recombination_rates= [1e-7, 1e-8, 1e-9]
v_sequence_length= 5e3
v_population_size= 1e4

def mrca_time_api(tree_sequence):
    ca_times = []
    for tree in tree_sequence.trees():
        ca_times.append(tree_sequence.tables.nodes[tree.root].time)
    return(min(ca_times))

def length_mrca(tree_sequence):
    min_ca_time= tree_sequence.max_root_time
    # min_interval= tree_sequence.sequence_length
    for tree in tree_sequence.trees():
        ca_time= tree_sequence.tables.nodes[tree.root].time
        interval= tree.interval.right-tree.interval.left
        if ca_time<= min_ca_time:
            min_ca_time= ca_time
            min_interval= interval
    return(min_interval)

def list_2_histogram(f_list, f_bin_number):
    bin_size= max(f_list)/f_bin_number
    bin_locations= [(idx+0.5)*bin_size for idx in range(f_bin_number)]
    bin_heights= [0 for idx in range(f_bin_number)]
    for entry in f_list:
        bin_idx= int(entry/bin_size)
        if bin_idx== f_bin_number:
            bin_idx= bin_idx-1
        bin_heights[bin_idx]= bin_heights[bin_idx]+1/len(f_list)/bin_size
    return [bin_locations, bin_heights]

def data_averaging(f_data, final_number):
    averaging_scope= max(f_data[0])/final_number
    final_iv1= [(idx+0.5)*averaging_scope for idx in range(final_number)]
    final_iv2= [0 for idx in range(final_number)]
    final_dv= [0 for idx in range(final_number)]
    points_collapsed= [0 for idx in range(final_number)]
    for idx in range(len(f_data[0])):
        point_idx= int(f_data[0][idx]/averaging_scope)
        if point_idx== final_number:
            point_idx= point_idx-1
        final_iv2[point_idx]= final_iv2[point_idx]+f_data[0][idx]
        final_dv[point_idx]= final_dv[point_idx]+f_data[1][idx]
        points_collapsed[point_idx]= points_collapsed[point_idx]+1
    for idx in range(final_number):
        if points_collapsed[idx]== 0:
            final_iv2[idx]= final_iv1[idx]
            final_dv[idx]= 0
        else:
            final_iv2[idx]= final_iv2[idx]/points_collapsed[idx]
            final_dv[idx]= final_dv[idx]/points_collapsed[idx]
    return [final_iv1, final_iv2, final_dv]

for v_recombination_rate in recombination_rates:
    times= []
    lengths= []
    for idx in range(trial_number):
        tree_sequence_sim = msprime.sim_ancestry(
            samples= sample_size,
            recombination_rate= v_recombination_rate,
            sequence_length= v_sequence_length,
            population_size= v_population_size,
            ploidy=1,
            record_full_arg=True)
        times.append(mrca_time_api(tree_sequence_sim))
        lengths.append(length_mrca(tree_sequence_sim))
    histogram= list_2_histogram(lengths, v_bin_number)
    with open("{:.1e}".format(v_recombination_rate)+"_histogram.txt", 'w') as txtfile:
        txtfile.write(str(histogram))
    averaged_data= data_averaging([times, lengths], v_bin_number)
    with open("{:.1e}".format(v_recombination_rate)+"_times.txt", 'w') as txtfile:
        txtfile.write(str(averaged_data))
