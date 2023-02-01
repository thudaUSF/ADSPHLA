import traceback
import pandas as pd
import numpy as np
from statsmodels.stats.proportion import proportions_ztest
import scipy.stats as stats
import rpy2.robjects as R
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
import csv
import os
import sys
import itertools
import math

"""
Overall, finds significant phenotypic differences of people w/ certain HLA within AD+ patients. Also compares HLA prevalence between people w/, w/o AD.
Pre-requisite files - ADHLApheno.csv - Have HLA and phenotype data with appropriate headers. 
    Ordered top to bottom by most HLA found, so more blanks near bottom, more likely removed, then ordered by name.
Workflow:
1. Create reference, filter to sampletype, Combine Braak 3,4.
2. Count frequency of HLA types. True frequency. (Removes duplicates if a person has the allele twice)
#4. Count frequency of V- or J-segments. True frequency.
#5. Appends HLA types to V- J-segments to create combination.
#6. Find frequency of combinations (not true frequency because one would have to drop specific combos that are duplicate in each case)
#8. Find survival log-rank test for HLA (Prints HLA significant cases) , VJ segment, and Combination. KM Plots for significant cases.
9. Print all significant cases. Outputs to excel file. #Prints arms for significant combinations.

"""

path = "/Users/thuda/Desktop/Research/14-ADHLA/2-Processing/" 

#Takes aliquot thing, ends up with filename, case, Submitter id (which is used to identify Origin) 
clinical = pd.read_csv(path + "ADHLApheno.csv",sep=',')
"""for column in clinical:
    if "HLA-" in column:
        clinical[column] = clinical[column].str.split(':').str[0]""" #removes everything after semicolon (e.g. no distinguishing between subtypes)
clinical.loc[(clinical['Braak'] < 3), 'Braak'] = 1
clinical.loc[(clinical['Braak'] == 4), 'Braak'] = 3
clinical = clinical.rename(columns={"SUBJID":"Case"})
sampletype = "All" #to change to All, comment the next line
#clinical = clinical[clinical["BODY_SITE"]==sampletype]
clinical = clinical.drop_duplicates(subset = "Case", keep = 'first', inplace = False)
clinical['APOE']= clinical['APOE'].map(str)
clinicalonly = clinical[["Case","SampleID","BODY_SITE",
                           'APOE',
                           "Age",'AD',"Braak"]] 
clinicalisnot = clinical[clinical['AD']== 0]
clinical = clinical[clinical['AD'] == 1]
HLA=clinical.copy()

#Count amount (frequency) of people with a certain HLA subtype. Get total count for all HLAs in both alleles - add total counts, subtract repetitions (people who had the same allele twice.)
#Count of HLA may be slightly different than a manual calculation because people can be typed differently depending on which duplicate removed (typing difference is usually between subtypes e.g. 3:01 --> 3:05)
frequentHLAtype = pd.DataFrame()
frequentHLAtypepairs = pd.DataFrame()
frequentHLA = pd.DataFrame()
frequentHLAtotal = pd.DataFrame(columns = ['HLA','frequency'])
HLAlistdf = pd.DataFrame()
for columnname in HLA:
    if "HLA-" in columnname and not "*" in columnname: #For every HLA type, don't do HLA-A'
        frequentHLAtypeSeries = HLA[columnname].value_counts() #Counts frequency of values in a list
        frequentHLAtype = pd.DataFrame({columnname:frequentHLAtypeSeries.index, 'count':frequentHLAtypeSeries.values}) #Convert to Dataframe; first column = HLA type, second column = frequency
        frequentHLAtypeSeries2 = HLA[columnname + "*"].value_counts()
        frequentHLAtype2 = pd.DataFrame({columnname:frequentHLAtypeSeries2.index, 'count*':frequentHLAtypeSeries2.values}) #Do above for second allele
        HLAtype = columnname.replace("HLA-","")
        #print(HLAtype)
        namepairs = HLAtype + 'pairs'
        HLA[namepairs] = HLA[columnname] + HLA[columnname + "*"] #Combine the alleles together
        frequentHLAtypepairsSeries = HLA[namepairs].value_counts() #Frequency of combined alleles
        frequentHLAtypepairs = pd.DataFrame({namepairs:frequentHLAtypepairsSeries.index, 'count pairs':frequentHLAtypepairsSeries.values}) #Convert to dataframe
        frequentHLAtypepairs[columnname] = frequentHLAtypepairs[namepairs].str.split(HLAtype, expand=False).str[1]
        frequentHLAtypepairs['2'] = frequentHLAtypepairs[namepairs].str.split(HLAtype, expand=False).str[2]
        #print(frequentHLAtypepairs.head())
        is_repeat = frequentHLAtypepairs[columnname] == frequentHLAtypepairs['2'] 
        frequentHLAtypepairs = frequentHLAtypepairs[is_repeat]
        frequentHLAtypepairs = frequentHLAtypepairs.reset_index(drop = True) #Take pairs from the total count where both alleles are the same only
        frequentHLAtypepairs[columnname] = HLAtype + frequentHLAtypepairs[columnname].astype(str) #Fixing string to make it match 
        frequentHLAtype = frequentHLAtype.merge(frequentHLAtype2, on = columnname,how='left') #Merge allele2
        frequentHLAtype = frequentHLAtype.merge(frequentHLAtypepairs[[columnname,'count pairs']], on = columnname,how='left') #Merge pairing numbers
        truefreqname = 'truefrequency HLA-' + HLAtype
        frequentHLAtype[truefreqname] = frequentHLAtype.fillna(0)['count'] + frequentHLAtype.fillna(0)['count*'] - frequentHLAtype.fillna(0)['count pairs']  #add the amounts between number of alleles, subtract pairs
        frequentHLAtype = frequentHLAtype.rename(columns={columnname:'HLA', truefreqname:'frequency'}) #match it so it can be added to a list of all HLA    
        frequentHLAtype = frequentHLAtype.sort_values(by='frequency',ascending=False)
        frequentHLAtype = frequentHLAtype.reset_index(drop = True)
        frequentHLAtotal = frequentHLAtotal.append(frequentHLAtype[['HLA','frequency']], ignore_index=True) #List of all HLAs in a list and real counts
    else:
        continue
is_boundary_HLA = frequentHLAtotal['frequency']>=40
ChecklistHLA = pd.DataFrame(frequentHLAtotal[is_boundary_HLA])

HLA = HLA.drop(columns=['SampleID','Apairs','Bpairs','Cpairs','DPB1pairs','DQB1pairs','DRB1pairs'])
conditions = []
allconditions = False
Cases = pd.DataFrame() #All Cases that fit the requirements
RCases = pd.DataFrame() #Replicative set of Cases, Randomly generates half from cases
Casesdrop = pd.DataFrame() #For Cases that fit requirements, remove duplicate people (for example: multiple reads of a HLA-TRB-V on same person that match) dropped
brainresults = []
ageresults = []
propresults = []
APOEresults = []
better=""

def filterdf(conditions,HLAwork):
    Used = HLAwork.copy()
    filters= Used['AD']!=2 #Baseline set to always true
    compfilters=filters
    VJ=False
    for condition in conditions:
        if "*" in condition:
            HLAtype = "HLA-" + condition.split('*')[0] #HLA-DPB
            column1 = HLAtype
            column2 = HLAtype + '*' 
            is_c1 = Used[column1]==condition
            is_c2 = Used[column2]==condition
            filters = filters&is_c1|is_c2 # is receptor found in Column 1 or 2 (HLA-A or HLA-A')
            compfilters=compfilters&~(is_c1|is_c2)
        if "notAPOE4" in condition:
            is_apo = Used['APOE'].str.contains("4")
            filters = filters&~is_apo
            compfilters = compfilters&~is_apo
        if "notAPOE2" in condition:
            is_apo = Used['APOE'].str.contains("2")
            filters = filters&~is_apo
            compfilters = compfilters&~is_apo

    Cases = Used[filters]

    #CompCases = Used[~Used['Case'].isin(Cases['Case'])] #Remove rows where the case is found in the Observed group, all variables with Comp --> Complementary Set
    CompCases = Used[compfilters&Used[HLAtype].notnull()]
    return Cases, CompCases

def run_MW(finallist,df,dfcomp):
    df["group"] = 1
    dfcomp["group"] = 0 
    frames = [df, dfcomp]
    pdf = pd.concat(frames)
    pdf = pdf[pdf.Braak != 0] #drop in the same way as KW
    pdf = pdf[["group","Braak"]]
    if df['Braak'].mean() > dfcomp['Braak'].mean():
        better="Group is worse"
    else:
        better="Group is better"
    with localconverter(R.default_converter + pandas2ri.converter):
        df = R.conversion.py2rpy(pdf)
    MW = R.r(r'''
        function(df) {
            options(warn=-1)
            df$Braak <- as.factor(df$Braak)
            levels(df$Braak)
            df$Braak <- ordered(df$Braak, levels = c("1","3","5","6"))

            # run mann-whitney for property
            pvalue <- wilcox.test(as.integer(Braak)~as.integer(group), data=df)$p.value
            }
    ''')
    rresults = MW(df)
    results=np.asarray(rresults)
    finallist.extend(results)
    finallist.append(better)
    return finallist

def run_ttest(finallist,df1,df2):
    results = stats.ttest_ind(df1['Age'], df2['Age'],nan_policy='omit') #perform calculation ignoring nan values
    finallist.extend([df1['Age'].mean(),df2['Age'].mean(),results.pvalue])
    if df1['Age'].mean() > df2['Age'].mean():
        better="Group is better"
    else:
        better="Group is worse"
    finallist.append(better)
    return finallist

def run_zproportion(finallist,df1,df2,receptor):
    casecount = np.array([len(df1.index), len(df1.index) + len(df2.index)])
    df3,df4 = filterdf([receptor],clinicalisnot)
    successcount = np.array([len(df1.index), len(df3.index)])
    totalcount = np.array([len(df1.index) + len(df2.index), len(df3.index) + len(df4.index)])
    
    proportioncase = successcount[0]/totalcount[0]
    proportioncontrol = successcount[1]/totalcount[1]
    stat, pval = proportions_ztest(successcount, totalcount)
    #print('{0:0.3f}'.format(pval))
    finallist.extend([proportioncase, len(df3.index),len(df4.index), proportioncontrol, pval])
    if proportioncase > proportioncontrol:
        better="Group is worse"
    else:
        better="Group is better"
    finallist.append(better)
    return finallist

###### GENERATE RESULTS OF ANALYSIS IN LIST FORMAT ##########    
def get_data(sampletype,conditions,df,dfcomp,test): 
    data = [sampletype]
    data.extend(conditions)
    finallist = [len(df.index),len(dfcomp.index)]
    if test == "MW":
        list = run_MW(finallist,df,dfcomp)
    if test == "ttest":
        finallist[0] = (df.Age.count())
        finallist[1] = (dfcomp.Age.count()) #since some don't have age
        list = run_ttest(finallist,df,dfcomp)
    if test == "zproportion":
        list = run_zproportion(finallist,df,dfcomp,conditions[0])
    data.extend(list)
    return data

#Testing for survival differences for all HLA in Checklist HLA
RHLA = HLA.sample(frac = .5) #Creates half replicative set that will be used for all run throughs of HLA (All R before variable name - Replicative Set)

for index, row in ChecklistHLA.iterrows():
    item = row[0] #item = HLA being checked from checklist e.g. A*02:01
    conditions = [item]
    df1,df2 = filterdf(conditions,HLA)

    df1_90per = df1.sample(frac = .5)
    df1_10per = df1[~df1['Case'].isin(df1_90per['Case'])]
    frames = [df2, df1_10per]
    df2_110per = pd.concat(frames)

    """if item == "DQB1*02:02":
        save = item.translate(str.maketrans('', '', '*:\'\/')) #Problem with saving file if there's asterisk/semicolon
        HLA.to_csv(path + "HLAwork.csv",index=False)
        df1.to_csv(path + save + sampletype + ".csv", index=False) #Save Cases with no duplicates
        df2.to_csv(path + save + sampletype + "comp.csv", index=False) #Save Cases with no duplicates    
        df2_110per.to_csv(path + save + sampletype + "df2_110per.csv", index=False) #Save Cases with no duplicates
        df1_10per.to_csv(path + save + sampletype + "df1_10per.csv", index=False) #Save Cases with no duplicates    
        df1_90per.to_csv(path + save + sampletype + "df1_90per.csv", index=False) #Save Cases with no duplicates"""

    data = get_data(sampletype,conditions,df1,df2,"MW")
    Rdata = get_data(sampletype,conditions,df1_90per,df2_110per,"MW")
    
    if data[4] < .05:
        brainresults.append(data)
        brainresults.append(Rdata)
        
Headings = ["Sampletype","Receptor","N - Group", "N - Outgroup","MW p-value", "Better"]

df_brainresults = pd.DataFrame(brainresults)
df_brainresults.columns = Headings
with pd.ExcelWriter("HLA_Braak_Age_Prop_APOE_" + sampletype + 'results.xlsx') as writer:
    df_brainresults.to_excel(writer, sheet_name='Braak (MW)')
