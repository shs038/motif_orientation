#when two motifs co-occur with each other, their orientation should skew from background 
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import scipy
from scipy import stats
from collections import Counter
%matplotlib inline
#count orientation
def count_orientation(motif_orientation):
    '''
    input:a pandas dataframe contains motifs orientation data
    output:a 3 rows pandas dataframe contains counts of each orientation of each motifs
    '''
    motifs = motif_orientation.columns.values #save motifs identity 
    #creaty a zeros matrix to sotre future orientation count data
    zero_data = np.zeros((3,motif_orientation.shape[1]),dtype=np.int)
    # conver zeros matrix to zeros dataframe
    count_frame = pd.DataFrame(zero_data, columns=motifs)
    count_frame.index=['+','-','?']
    for i in range (motif_orientation.shape[1]):#loop to count the orientations of all motifs 
        one_motif=motif_orientation.ix[:,i].values#retrieve orientations of one motif
        one_motif=list(one_motif)#convert to list
        c=Counter(one_motif)#count each orientation
        c=dict(c)#convert to dictionary
        c=pd.DataFrame.from_dict(c,orient='index')#convert to dataframe
        count_frame.ix[:,i]=c.ix[:,0]#store orientation count in zeros dataframe 
    return count_frame
#read in orientation data
motifs_orientation_df = pd.read_excel('/home/shs038/chr1_motifs/chr1_orientation.xlsx',
                                      parse_cols="B:GO")#load part of file contain data
motifs_orientation_count=count_orientation(motifs_orientation_df)
#transpose the count dataframe 
motifs_orientation_count = motifs_orientation_count.T
#make stacked bar chart 
motifs_orientation_count.plot.barh(stacked=True,figsize=(24,60),)
#plot each motif's orientation
def motif_orientation_pie(motif):
    '''
    input: a string of motif identity
    output: a pie chart of the occurence of different orientations for the motif
    '''
    motifs_to_plot=motifs_orientation_count.ix[motif,0:2]#plot only '+' and '-'
    motifs_to_plot.plot.pie(autopct='%.2f', fontsize=10,figsize=(6,6))
motifs=motifs_orientation_df.columns.values
#count how many times two motifs that co-occur with each other both have sense orientation 
def count_both_sense(motif_orientation):
    '''
    input:a pandas dataframe contains motifs orientation data
    output:a pandas dataframe contains how many times each pair of motifs both have 
            + orientation.
    '''
    motifs = motif_orientation.columns.values #save motifs identity 
    #creaty a zeros matrix to sotre future orientation count data
    zero_data = np.zeros((motif_orientation.shape[1],motif_orientation.shape[1]),dtype=np.int)
    # conver zeros matrix to zeros dataframe
    count_frame = pd.DataFrame(zero_data, columns=motifs)
    count_frame.index=motifs
    for i in range (motif_orientation.shape[1]-1):
        #find the loci where the motif occur with sense orientation
        logical_col_i=motif_orientation.ix[:,i]=='+' 
        for j in range (i+1,motif_orientation.shape[1]):
            #find the loci where the motif occur with sense orientation
            logical_col_j=motif_orientation.ix[:,j]=='+'
            #find the loci where both of the motifs occur with sense orientation
            logical_input=1*logical_col_i+1*logical_col_j==2
            #count how many times both of the motifs occur with sense orientation
            input=np.sum(1*logical_input)
            count_frame.ix[i,j]=count_frame.ix[i,j]+input
    return count_frame
both_sense=count_both_sense(motifs_orientation_df)
#count how many times two motifs that co-occur with each other both have antisense orientation 
def count_both_antisense(motif_orientation):
    '''
    input:a pandas dataframe contains motifs orientation data
    output:a pandas dataframe contains how many times each pair of motifs both have 
            - orientation.
    '''
    motifs = motif_orientation.columns.values #save motifs identity 
    #create a zeros matrix to sotre future orientation count data
    zero_data = np.zeros((motif_orientation.shape[1],motif_orientation.shape[1]),dtype=np.int)
    # conver zeros matrix to zeros dataframe
    count_frame = pd.DataFrame(zero_data, columns=motifs)
    count_frame.index=motifs
    for i in range (motif_orientation.shape[1]-1):
        logical_col_i=motif_orientation.ix[:,i]=='-'
        for j in range (i+1,motif_orientation.shape[1]):
            logical_col_j=motif_orientation.ix[:,j]=='-'
            logical_input=1*logical_col_i+1*logical_col_j==2
            input=np.sum(1*logical_input)
            count_frame.ix[i,j]=count_frame.ix[i,j]+input
    return count_frame
both_antisense=count_both_antisense(motifs_orientation_df)
#count how many times two motifs co-occur with each other have sense/antisense orientation 
def count_sense_antisense(motif_orientation):
    '''
    input:a pandas dataframe contains motifs orientation data
    output:a pandas dataframe contains how many times each pair of motifs  have 
            +/- orientation.
    '''
    motifs = motif_orientation.columns.values #save motifs identity 
    #create a zeros matrix to sotre future orientation count data
    zero_data = np.zeros((motif_orientation.shape[1],motif_orientation.shape[1]),dtype=np.int)
    # conver zeros matrix to zeros dataframe
    count_frame = pd.DataFrame(zero_data, columns=motifs)
    count_frame.index=motifs
    for i in range (motif_orientation.shape[1]-1):
        logical_col_i=motif_orientation.ix[:,i]=='+'
        for j in range (i+1,motif_orientation.shape[1]):
            logical_col_j=motif_orientation.ix[:,j]=='-'
            logical_input=1*logical_col_i+1*logical_col_j==2
            input=np.sum(1*logical_input)
            count_frame.ix[i,j]=count_frame.ix[i,j]+input
    return count_frame
sense_antisense=count_sense_antisense(motifs_orientation_df)
#count how many times two motifs co-occur with each other have antisense/sense orientation 
def count_antisense_sense(motif_orientation):
    '''
    input:a pandas dataframe contains motifs orientation data
    output:a pandas dataframe contains how many times each pair of motifs  have 
            -/+ orientation.
    '''
    motifs = motif_orientation.columns.values #save motifs identity 
    #creaty a zeros matrix to sotre future orientation count data
    zero_data = np.zeros((motif_orientation.shape[1],motif_orientation.shape[1]),dtype=np.int)
    # conver zeros matrix to zeros dataframe
    count_frame = pd.DataFrame(zero_data, columns=motifs)
    count_frame.index=motifs
    for i in range (motif_orientation.shape[1]-1):
        logical_col_i=motif_orientation.ix[:,i]=='-'
        for j in range (i+1,motif_orientation.shape[1]):
            logical_col_j=motif_orientation.ix[:,j]=='+'
            logical_input=1*logical_col_i+1*logical_col_j==2
            input=np.sum(1*logical_input)
            count_frame.ix[i,j]=count_frame.ix[i,j]+input
    return count_frame
antisense_sense=count_antisense_sense(motifs_orientation_df)
def reshape_both_sense(co_orientation):
    '''
    input:a 196*196 pandas dataframe contains how many times each pair of motifs both have 
            + orientation.
    output: a reshaped pandas dataframe contains the same data but with 
            columns: orientation and counts, index:motif pair
    '''
    coorientation_dic={}#create an empty dictionary
    #loop in part of count data that contain meaning counting
    for i in range (co_orientation.shape[1]-1):
        for j in range (i+1,co_orientation.shape[1]):
            motif1 = motifs[i]
            motif2 = motifs[j]
            #put motif pair and count into the empty dictinory
            coorientation_dic[motif1,motif2]=co_orientation.ix[i,j]
    #convert the dictionary to dataframe
    reshaped_frame = pd.DataFrame(coorientation_dic.items(), 
                                  columns=['pair','Count'],
                                  index=dict.keys(coorientation_dic))
    reshaped_frame['Orientation']='+/+' #add orientation to the dataframe
    del reshaped_frame['pair']#remove pair as it is the same as index
    reshaped_frame=reshaped_frame[['Orientation','Count']]# put Orientation in front of Count
    return reshaped_frame
both_sense_reshaped=reshape_both_sense(both_sense)
def reshape_both_antisense(co_orientation):
    '''
    input:a 196*196 pandas dataframe contains how many times each pair of motifs both have 
            - orientation.
    output: a reshaped pandas dataframe contains the same data but with 
            columns: orientation and counts, index:motif pair
    '''
    coorientation_dic={}
    for i in range (co_orientation.shape[1]-1):
        for j in range (i+1,co_orientation.shape[1]):
            motif1 = motifs[i]
            motif2 = motifs[j]
            coorientation_dic[motif1,motif2]=co_orientation.ix[i,j]
    reshaped_frame = pd.DataFrame(coorientation_dic.items(), 
                                  columns=['pair','Count'],
                                  index=dict.keys(coorientation_dic))
    reshaped_frame['Orientation']='-/-'
    del reshaped_frame['pair']
    reshaped_frame=reshaped_frame[['Orientation','Count']]
    return reshaped_frame
both_antisense_reshaped=reshape_both_antisense(both_antisense)
def reshape_sense_antisense(co_orientation):
    '''
    input:a 196*196 pandas dataframe contains how many times each pair of motifs have 
            +/- orientation.
    output: a reshaped pandas dataframe contains the same data but with 
            columns: orientation and counts, index:motif pair
    '''
    coorientation_dic={}
    for i in range (co_orientation.shape[1]-1):
        for j in range (i+1,co_orientation.shape[1]):
            motif1 = motifs[i]
            motif2 = motifs[j]
            coorientation_dic[motif1,motif2]=co_orientation.ix[i,j]
    reshaped_frame = pd.DataFrame(coorientation_dic.items(), 
                                  columns=['pair','Count'],
                                  index=dict.keys(coorientation_dic))
    reshaped_frame['Orientation']='+/-'
    del reshaped_frame['pair']
    reshaped_frame=reshaped_frame[['Orientation','Count']]
    return reshaped_frame
sense_antisense_reshaped=reshape_sense_antisense(sense_antisense)
def reshape_antisense_sense(co_orientation):
    '''
    input:a 196*196 pandas dataframe contains how many times each pair of motifs have 
            -/+ orientation.
    output: a reshaped pandas dataframe contains the same data but with 
            columns: orientation and counts, index:motif pair
    '''
    coorientation_dic={}
    for i in range (co_orientation.shape[1]-1):
        for j in range (i+1,co_orientation.shape[1]):
            motif1 = motifs[i]
            motif2 = motifs[j]
            coorientation_dic[motif1,motif2]=co_orientation.ix[i,j]
    reshaped_frame = pd.DataFrame(coorientation_dic.items(), 
                                  columns=['pair','Count'],
                                  index=dict.keys(coorientation_dic))
    reshaped_frame['Orientation']='-/+'
    del reshaped_frame['pair']
    reshaped_frame=reshaped_frame[['Orientation','Count']]
    return reshaped_frame
antisense_sense_reshaped=reshape_antisense_sense(antisense_sense)
#concatenate data for all four orientation
frames = [both_sense_reshaped,both_antisense_reshaped,
          sense_antisense_reshaped,antisense_sense_reshaped]
co_occur_orientation = pd.concat(frames)
#plot without Count=0
antisense_sense_reshaped_nonz=antisense_sense_reshaped[antisense_sense_reshaped['Count']>0]
sense_antisense_reshaped_nonz=sense_antisense_reshaped[sense_antisense_reshaped['Count']>0]
both_sense_reshaped_nonz=both_sense_reshaped[both_sense_reshaped['Count']>0]
both_antisense_reshaped_nonz=both_antisense_reshaped[both_antisense_reshaped['Count']>0]
sns.distplot(antisense_sense_reshaped_nonz['Count'])
sns.distplot(sense_antisense_reshaped_nonz['Count'])
sns.distplot(both_sense_reshaped_nonz['Count'])
sns.distplot(both_antisense_reshaped_nonz['Count'])
#normalize orientation count
fd_add=both_sense_reshaped['Count'].add(both_antisense_reshaped['Count'], fill_value=0)
fd_add=fd_add.add(sense_antisense_reshaped['Count'], fill_value=0)
fd_add=fd_add.add(antisense_sense_reshaped['Count'], fill_value=0)
co_occur_orientation['Total']=fd_add
co_occur_orientation = co_occur_orientation[co_occur_orientation.Total != 0]
co_occur_orientation['Normalized_Count']=co_occur_orientation['Count']/co_occur_orientation['Total']
#subset of co_occur_orientation for each orientation pair without Normalized_Count=0
co_occur_sense=co_occur_orientation[co_occur_orientation['Orientation']=='+/+']
co_occur_sense_nz=co_occur_sense[co_occur_sense['Normalized_Count']>0]
co_occur_antisense=co_occur_orientation[co_occur_orientation['Orientation']=='-/-']
co_occur_antisense_nz=co_occur_antisense[co_occur_antisense['Normalized_Count']>0]
co_occur_sense_antisense=co_occur_orientation[co_occur_orientation['Orientation']=='+/-']
co_occur_sense_antisense_nz=co_occur_sense_antisense[co_occur_sense_antisense['Normalized_Count']>0]
co_occur_antisense_sense=co_occur_orientation[co_occur_orientation['Orientation']=='-/+']
co_occur_antisense_sense_nz=co_occur_antisense_sense[co_occur_antisense_sense['Normalized_Count']>0]
sns.distplot(co_occur_sense_nz['Normalized_Count'])
plt.show()
sns.distplot(co_occur_antisense_nz['Normalized_Count'])
plt.show()
sns.distplot(co_occur_sense_antisense_nz['Normalized_Count'])
plt.show()
sns.distplot(co_occur_antisense_sense_nz['Normalized_Count'])
# Given a set of loci L where motif J is in a given orientation O, does that subset of L where motf Z is present have a bias towards +/-
# Null hypothesis: orientation of motif J and that of motif Z are independent at loci L.
def Fisher_test(sense_sense_df,antisense_antisense_df,sense_antisense_df,antisense_sense_df):
    # create a  2*2 dataframe for fisher exact test
    Fisher_df=pd.DataFrame(data=np.zeros((2,2),dtype=np.int),
                       columns=['motif_1_sense','motif_1_antisense'], 
                       index=['motif_2_sense','motif_2_antisense'])
    motifpair=sense_sense_df.index.values
    P_dic={}#create a dictionary to store p value
    for pair in motifpair:
        Fisher_df['motif_1_sense']['motif_2_sense']=sense_sense_df.loc[(pair,),'Count']
        Fisher_df['motif_1_antisense']['motif_2_sense']=sense_antisense_df.loc[(pair,),'Count']
        Fisher_df['motif_1_sense']['motif_2_antisense']=antisense_sense_df.loc[(pair,),'Count']
        Fisher_df['motif_1_antisense']['motif_2_antisense']=antisense_antisense_df.loc[(pair,),
                                                                                       'Count']
        oddsratio,p=stats.fisher_exact(Fisher_df)
        P_dic[pair]=p
    P_df=pd.DataFrame.from_dict(P_dic,orient='index')#convert p dictionary to a dataframe
    return P_df
P_dataframe=Fisher_test(co_occur_sense,co_occur_antisense,
                        co_occur_sense_antisense,co_occur_antisense_sense)
sns.distplot(P_dataframe)
#find the motif pair whose p <0.05, corret by 195
P_correct_by_times=P_dataframe[P_dataframe.ix[:,0]<(0.05/195)]
sns.distplot(P_correct_by_times)
significant_pairs = P_correct_by_times.index.values
for sp in significant_pairs:
    if sp[0] == 'rel' or sp[1] == 'rel':
        print(sp)
co_occur_pairs = pd.read_csv('/home/shs038/chr1_motifs/sigpair.tsv', sep='\t')
del co_occur_pairs['Unnamed: 0']
Motifpair = P_correct_by_times.index.values
for mp in Motifpair:
    for j in range (len(co_occur_pairs)):
        if mp[0] == co_occur_pairs.loc[j,'motif1'] and mp[1] == co_occur_pairs.loc[j,'motif2']:
            P_correct_by_times.loc[(mp,),:]=100
        elif mp[0] == co_occur_pairs.loc[j,'motif2'] and mp[1] == co_occur_pairs.loc[j,'motif1']:
            P_correct_by_times.loc[(mp,),:]=100
P_correct_by_times=P_correct_by_times[P_correct_by_times.ix[:,0]!=100]
