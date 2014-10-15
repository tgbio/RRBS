#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 14 11:06:54 2014

@author: tg
"""
#%%
from Bio import Restriction 
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from BCBio import GFF
from Bio.Alphabet import generic_dna
#from BCBio.GFF import GFFExaminer
import csv
from sys import exit
#import pdb
#%%
#################################################
#Define relevant enzymes and their cut sites
#################################################
print "Defining relevant enzymes"
#Reference: http://biopython.org/DIST/docs/cookbook/Restriction.html

rb = Restriction.RestrictionBatch([])

indiv_enzymes = ["BsaWI",  "BssSI", "BsoBI"]
double_digs = [("BsaWI",  "BssSI"), ("BsaWI", "BsoBI"), ("BssSI", "BsoBI")]
single_double_digs = indiv_enzymes + double_digs
for c in indiv_enzymes:
    assert c in Restriction.CommOnly
    rb.add(c)

#%%
########################################
#Determine the number of bins per contig
########################################
print "Determining the number of bins per contig"
binsize = 10000

print "Binsize =", binsize

print "Defining relevant enzymes and their cut sites for each contig"

c_filelist = ["Lsat.1.v4.ScaffoldsSequence.All.fasta_biggerthan_1000.5.fasta"] #N.5 is TEMP for testing, note that the contigs must be at least 1000 bp and the bigest contig is 3.1 million bp

contig_bin_dict = {} #a dictionary with keys = contig ids and values = 1 if that contig's sequence is 10k bp or less, or else length of the sequnce/10k, which will later become the number of bins
#temp_contig_limit = 0
for f in c_filelist: #for contig file in c_filelist, normally will be just one file with all contigs
    for seq_record in SeqIO.parse(f, "fasta"): # for each contig (fasta record) in contig file
        
        lencontig = len(str(seq_record.seq)) # add the contig name (fasta header) to contig_id_dict, a dictionary with key = contig name, value = length of sequence
        if lencontig <= binsize: # if the total length of the contig is the same size or smaller than the binsize, its not worth binning it
            contig_bin_dict[seq_record.id] = 1
        else:
            contig_bin_dict[seq_record.id] = lencontig/binsize
#%%
#################################################
#Define interesting regions of the lettuce genome
#################################################
print "Begining section to define interesting regions of the lettuce genome."
#Reference:http://biopython.org/wiki/GFF_Parsing

#Interesting regions defined by these gff files:
#Lsat.1.v4.geneticMarkers.Lsativa-Lserriola_core_map.gff - GFF file with the genetic markers present in the L. sativa- L. serriola; gff_id = all contig names, 'gff_source_type': {('Salinas_x_Serriola_2011_Good', 'Map_Locations')
#Lsat.1.v4.alignment.allPlantProteins.gff - GFF exonerate alignment of plant protein database from NCBI; gff_id = all the contigs names, gff_source_type ('protein2genome', 'protein_match'):164235
#Lsat.1.v4.repeatElements.gff - GFF file with the predicted repeated elements
#Lsat.1.v4.alignment.lettuceNCBI.gff - GFF exonerate alignment of lettuce mRNA sequences at NCBI; gff_source = est2genome, gff_type = "expressed_sequence_match"
#Lsat.1.v4.geneModels.gff - GFF file with the predicted gene models

#limit_info is a dictionary defining the interesting fields from each gff 
#limit_info = {"Lsat.1.v4.geneticMarkers.Lsativa-Lserriola_core_map.gff": {'gff_source_type':[('Salinas_x_Serriola_2011_Good', 'Map_Locations')]}}#, "Lsat.1.v4.alignment.lettuceNCBI.gff":{"gff_source_type":[("est2genome", "expressed_sequence_match")]}, "Lsat.1.v4.alignment.allPlantProteins.gff":{"gff_source_type":[('protein2genome', 'protein_match')]}, "Lsat.1.v4.repeatElements.gff": {'gff_source_type': [('RepeatMasker_v4', 'Transposon'),('RepeatProteinMask_V3', 'TEprotein'),('TRF_v4', 'TandemRepeat')]}}
#limit_info = {"Lsat.1.v4.repeatElements.gff": {'gff_source_type': [('RepeatMasker_v4', 'Transposon'),('RepeatProteinMask_V3', 'TEprotein'),('TRF_v4', 'TandemRepeat')]}}
limit_info = {"Lsat.1.v4.alignment.lettuceNCBI.gff":{"gff_source_type":[("est2genome", "expressed_sequence_match")]}, "Lsat.1.v4.repeatElements.gff": {'gff_source_type': [('RepeatMasker_v4', 'Transposon'),('RepeatProteinMask_v4', 'TEprotein'),('TRF_v4', 'TandemRepeat')]}}
limit_info = {"Lsat.1.v4.alignment.lettuceNCBI.gff":{"gff_source_type":[("est2genome", "expressed_sequence_match")]}}
all_contig_dict = {x:{z:"contig not in gff" for z in limit_info} for x in contig_bin_dict} # a dictionary with keys = contig ids and values = a dictionary with keys = gff filename and value = a list of as many lists as there are bins for that contig, with each of the bin lists containing lists of [begining position of feature, ending position of feature]


for in_file in limit_info: # loop over each .gff file
    print in_file
    in_handle = open(in_file) # open up the .gff file
    print "limit_info[in_file]", limit_info[in_file] # print the limits (the interesting fields) for that .gff
    temp = 0 #TEMP for limiting the number of records to analyze
    for contig in GFF.parse(in_handle, limit_info=limit_info[in_file]): #for each contig in the gff file having a limit value
        if contig.id in all_contig_dict:   
            all_contig_dict[contig.id][in_file] = [[]]*(contig_bin_dict[contig.id]+1) #a list of as many lists as there are bins for that contig, with each of the bin lists containing lists of [begining position of feature, ending position of feature]
            temp_feat_list = [[]]*(contig_bin_dict[contig.id]+1) #a list of as many lists as there are bins for that contig, with each of the bin lists containing lists of [begining position of feature, ending position of feature]
            if contig.id == "Lsat_1_v4_g_0_1178": #TEMP
                print "****FIRST len(all_contig_dict[contig.id][in_file])", len(all_contig_dict[contig.id][in_file])#TEMP
        
        else:
            print "GFF contig.id", contig.id, "not in contig fasta"
            continue
        print "Number of features on contig", len(contig.features)        
        #print "First feature", contig.features[0]

        #this section checks for overlapping positions, expands boundaries of feature to include overlaps if necessary, or adds non-overlapping features
        for feat in contig.features: #for each feature 
            binassign = feat.location.start/binsize            
            minf = min(int(feat.location.start), int(feat.location.end))
            maxf = max(int(feat.location.start), int(feat.location.end))
            if contig.id in all_contig_dict: # if the contig.id is in the all_contig_dict
                if len(temp_feat_list[binassign]) == 0:
                    temp_feat_list[binassign].append([minf, maxf])
                else:
                    for x in temp_feat_list[binassign]:
                        if minf in range(x[0], x[1]) and maxf in range(x[0], x[1]):
                            break
                        elif minf in range(x[0], x[1]) and maxf > x[1]:
                            temp_feat_list[binassign].remove(x)
                            temp_feat_list[binassign].append([x[0], maxf])
                            break
                        elif maxf in range(x[0], x[1]) and minf < x[0]:
                            temp_feat_list[binassign].remove(x)
                            temp_feat_list[binassign].append([minf, x[1]])
                            break
                        elif minf < x[0] and maxf > x[1]:
                            temp_feat_list[binassign].remove(x)
                            temp_feat_list[binassign].append([minf, maxf])
                            break
                        else:
                            all_contig_dict[contig.id][in_file][binassign].append([minf, maxf])
                all_contig_dict[contig.id][in_file] = temp_feat_list
        if temp >= 50:
            break
        temp +=1

    in_handle.close() # close the .gff file
    
print all_contig_dict
#exit()
    
#%%
############################################################################
#Categorize useful fragments by their location relative to genome features
############################################################################
    
def makesortfragments(cutter, all_sites, featuredigest, targetdigseqs, binsize, seq_record):

    #################################
    #Find, count the number and length and C content of the interior re fragments of the contig
    #################################
    #print len(all_sites)
    approx_fragments = [seq_record.seq[all_sites[i]:all_sites[i+1]] for i in range(len(all_sites)-1)] #approximate fragment as not taking into account any nt appearing before the cut site
    #print "approx_fragments", approx_fragments    
    
    targetdigseqs[cutter] = targetdigseqs[cutter] + [i for i in approx_fragments if len(i) in range(200, 500)]
    #print len(targetdigseqs[cutter])
    #break #TEMP
            
            
    interior_fragments = approx_fragments #the first and last fragments exclueded as they are not biologically meaningful; the start and end position are artifacts of the read assembly not ends of biologically occuring molecules
    #print "interior_fragments", interior_fragments
    #print 'len(interior_fragments)', len(interior_fragments)
    #print "len(all_sites)-1", len(all_sites)-1
            
        
    assert len(interior_fragments) == len(all_sites)-1 #the number of interior fragments should be equal to the number of cut sites minus one, as the total number of fragments is equal to the number of cut sites plus one
            
    fragfacts = [[m.count("C"), len(m)] for m in interior_fragments] # as list of lists, where each internal list contains in the first position a count of the number of C's in that internal fragment, and the second position contains the length of that internal fragment

    #print "range(len(cutsites)-1)", range(len(cutsites)-1)
    #if type(cutter) == str:
    #    print all_sites
    #    print approx_fragments
        
    for r in range(len(all_sites)-1): #where len(cutsites)-1 == len(interior_fragments)
        #print fragfacts[r]                
        bpos = int(all_sites[r]) 
        epos = int(all_sites[r+1])
        #print bpos, epos                
        
        fragfacts[r].append([bpos, epos]) #elements of fragfacts = [count C in frag, len frag, [5' cut site, 3' cut site]]
        #print fragfacts[r]                
        #break    
        if bpos/binsize == epos/binsize:
            lookupbin = [bpos/binsize]
        else:
            lookupbin = range(min(bpos/binsize, epos/binsize), max(bpos/binsize, epos/binsize)+1)
                
        for in_file in limit_info:
            if all_contig_dict[seq_record.id][in_file] == "contig not in gff":
                #print "type(cutter)", type(cutter)
                #print "contig not in gff"
                break
            print in_file
            for l in lookupbin:
                print type(l)
                print l, "bin"
                print all_contig_dict[seq_record.id][in_file][l]
                for region in all_contig_dict[seq_record.id][in_file][l]: ## a dictionary with keys = contig ids and values = a dictionary with keys = gff filename and value = a list of as many lists as there are bins for that contig, with each of the bin lists containing lists of [begining position of feature, ending position of feature]
                    print "#####################all_contig_dict[seq_record.id][in_file][l]",all_contig_dict[seq_record.id][in_file][l]                  
                    print "#####################region in all_contig_dict[seq_record.id][in_file][l]", region
                    #exit()
                    print region                   
                    if bpos in range(region[0], region[1]+1) or epos in range(region[0], region[1]+1) or min(bpos, epos) < min(region[0], region[1]) and max(bpos, epos) < max(region[0], region[1]):
                        if type(cutter) == str:#test, all the single digest fragments seem to be falling mostly within this first bin - is that realistic?
                            print cutter
                            print "begining cut site", bpos
                            print "ending cut site", epos
                            print "feature begin, end", region[0], region[1]
                        
                        featuredigest[in_file][cutter][0].append(fragfacts[r]) #if its in the feature region, append elements of fragfacts for that cutsite = [count C in frag, len frag, [5' cut site, 3' cut site]] to the first list
                        print "cutsite within feature"
                        break
                else: #excuted if for loop over all regions executes without breaking (i.e. cut site was not within any feature region)
                    for region in all_contig_dict[seq_record.id][in_file][l]:
                        if bpos in range(region[0]-1000, region[1]+1001) or epos in range(region[0]-1000, region[1]+1001):
                            featuredigest[in_file][cutter][1].append(fragfacts[r]) #if its not in the feature region, but within 1kb +/- append elements of fragfacts for that cutsite = [count C in frag, len frag, [5' cut site, 3' cut site]] to the second list
                            print "cutsite not in feature, within 1 kb +-"                                    
                            break
                    else:  #excuted if for loop over all regions executes without breaking (i.e. cut site was not within +/- 1kb of any feature region)
                        for region in all_contig_dict[seq_record.id][in_file][l]:                            
                            if bpos in range(region[0]-5000, region[1]+5001) or epos in range(region[0]-5000, region[1]+5001):
                                featuredigest[in_file][cutter][2].append(fragfacts[r]) #if its not in the feature region, but within 5kb +/- append elements of fragfacts for that cutsite = [count C in frag, len frag, [5' cut site, 3' cut site]] to the third list
                                print "cutsite not in feature, within 5 kb +-"                                 
                                break
                        else: #excuted if for loop over all regions executes without breaking (i.e. cut site was not within +/- 5kb of any feature region)
                            for region in all_contig_dict[seq_record.id][in_file][l]:
                                if bpos in range(region[0]-10000, region[1]+10001) or epos in range(region[0]-10000, region[1]+10001):
                                    featuredigest[in_file][cutter][3].append(fragfacts[r]) #if its not in the feature region, but within 10kb +/- append elements of fragfacts for that cutsite = [count C in frag, len frag, [5' cut site, 3' cut site]] to the fourth list
                                    print "cutsite not in feature, within 10 kb +-" 
                                    break
                            else: #excuted if for loop over all regions executes without breaking (i.e. cut site was not within +/- 10kb of any feature region)
                                print "cutsite not in feature or within 10 kb +-"                                 
                                featuredigest[in_file][cutter][4].append(fragfacts[r]) #if its not in the feature region nor within 10kb of the region append to the fifth list

    return featuredigest, targetdigseqs
    
print "Begining section to categorize useful fragments by their location relative to genome features."

featuredigest = {a:{x:[[], [], [], [], []] for x in single_double_digs} for a in limit_info.keys()}
 # a dictionary with keys = gff filenames, values = a dictionary with keys = tuple of double digest enzymes, values = a list of five lists of the features [count C in frag, len frag, [5' cut site, 3' cut site]] of each fragment, with one list for each proximity to the relevant gff feature (in feature, within 1kb, within 5kb, within 10kb, greater than 10kb)
targetdigseqs = {i:[] for i in single_double_digs} # a dictionary with keys = tuple of double digest enzymes and values = digestion fragments between 200 and 500 bp in length

#c_filelist = ["Lsat.1.v4.ScaffoldsSequence.All.fasta_biggerthan_1000.fasta"] 
c_filelist = ["Lsat.1.v4.ScaffoldsSequence.All.fasta_biggerthan_1000.5.fasta"] #temp for testing
#c_filelist = ["BsaWIBssSIBsoBI_double_digest_test2.fasta"]

for f in c_filelist: # for contig fasta file
    #record_counter = 0 #TEMP for testing
    for seq_record in SeqIO.parse(f, "fasta"): #get the contig sequence
        
        #record_counter +=1 #TEMP
        #if record_counter >=5000:  #TEMP
        #    break #TEMP
        #print "seq_record.id", seq_record.id
        #seqname = seq_record.id
        ana = Restriction.Analysis(rb, seq_record.seq, linear=True) #RestrictionBatch can give you a dictionary with the sites for all the enzymes in a RestrictionBatch
        #print "ana.with_sites()", ana.with_sites()
        #print "len(ana.with_sites())", len(ana.with_sites())
        if len(ana.with_sites()) == 1: #if there is only one enzyme that cuts within a fragement
            #print "len(ana.with_sites()) == 1"
            for x in ana.with_sites(): #loop over that one key in dictionary
               # print "x", x
                #print "len(ana.with_sites()[x])", len(ana.with_sites()[x])
                if len(ana.with_sites()[x]) > 1:#ask if the enzyme has more than one cut site on the contig
                    all_sites = ana.with_sites()[x]
                    all_sites.sort()
                    print "single all_sites", all_sites
                    cutter = str(x) 
                    #print "cutter, all_sites", cutter, all_sites
                    featuredigest, targetdigseqs = makesortfragments(cutter, all_sites, featuredigest, targetdigseqs, binsize, seq_record )
                    #featuredigest, targetdigseqs = makesortfragments(cutter, all_sites, seq_record.id)#make a call to a function for the fragment generation, calculation section
                    
        else: #there are more than 2 enzymes with cutsites on the contig
            enzyme_wrangler = {str(re):ana.with_sites()[re] for re in ana.with_sites()}
            #print "enzyme_wrangler", enzyme_wrangler
            for pair in double_digs:
                #print "pair[0]", pair[0]
                #print "pair[1]", pair[1]
                if pair[0] in enzyme_wrangler and pair[1] in enzyme_wrangler:                    
                    cutter = pair
                    #print "cutter", cutter
                    all_sites = enzyme_wrangler[pair[0]] + enzyme_wrangler[pair[1]]
                    #print all_sites
                    all_sites = [i for i in set(all_sites)]
                    #print "double all_sites", all_sites
                    all_sites.sort()
                    
                    if len(all_sites) < 2: #if there are not at least 2 cut sites on the fragment, no interior fragments will be generated, so continue on with the
                        continue
                    #print "cutter, all_sites", cutter, all_sites
                    featuredigest, targetdigseqs = makesortfragments(cutter, all_sites, featuredigest, targetdigseqs, binsize, seq_record )
                    #featuredigest, targetdigseqs = makesortfragments(cutter, all_sites, seq_record.id)#make a call to a function for the fragment generation, calculation section


################################################################################
##Consolidate singletons with appropriate double digests in featuredigest and targetdigseqs
################################################################################
#TODO something is still wrong with the featuredigst
doublefeaturedigest = {a:{} for a in limit_info.keys()} #create a new dictionary, with key = double digest tuple, value = a list 
for gff in featuredigest: ## a dictionary with keys = gff filenames, values = a dictionary with keys = tuple of double digest enzymes, values = a list of five lists of the features [count C in frag, len frag, [5' cut site, 3' cut site]] of each fragment, with one list for each proximity to the relevant gff feature (in feature, within 1kb, within 5kb, within 10kb, greater than 10kb)
    for enzyme_entry in featuredigest[gff]:
        print "things about featuredigest"        
        print enzyme_entry, len(featuredigest[gff][enzyme_entry])
        for r in featuredigest[gff][enzyme_entry]:
            print len(r)
        print "\n"
        if type(enzyme_entry) is tuple: #if the entry is a tuple (double digest)
            if enzyme_entry[0] in featuredigest[gff]: #if the first enzyme is in feature digest (if a contig had only a single enzyme's cut sites)
                assert len(featuredigest[gff][enzyme_entry]) == len(featuredigest[gff][enzyme_entry[0]])              
                for i in range(len(featuredigest[gff][enzyme_entry[0]])):
                    featuredigest[gff][enzyme_entry][i] = featuredigest[gff][enzyme_entry][i] + featuredigest[gff][enzyme_entry[0]][i] #then add the single digest fragments to the double digest fragments
            if enzyme_entry[1] in featuredigest[gff]: #if the second enzyme is in feature digest (if a contig had only a single enzyme's cut sites)
                assert len(featuredigest[gff][enzyme_entry]) == len(featuredigest[gff][enzyme_entry[1]])              
                for i in range(len(featuredigest[gff][enzyme_entry[1]])):
                    featuredigest[gff][enzyme_entry][i] = featuredigest[gff][enzyme_entry][i] + featuredigest[gff][enzyme_entry[1]][i] #then add the single digest fragments to the double digest fragments
            doublefeaturedigest[gff][enzyme_entry] = featuredigest[gff][enzyme_entry]
            
print "doublefeaturedigest"
for d in doublefeaturedigest[gff]:
    print d, len(doublefeaturedigest[gff][d]), type(d)
    for q in doublefeaturedigest[gff][d]:
        print len(q)
    print '\n'

print "Begin double target section"
# Create a dictionary to contain both the single and double digest fragments for the enzyme pairs
doubletargetdigseqs = {} # a dictionary with keys = tuple of double digest enzymes or single enzyme and values = single and double digestion fragments between 200 and 500 bp in length
for enzyme_entry in targetdigseqs: # a dictionary with keys = tuple of double digest enzymes or single enzyme and values = double digestion fragments between 200 and 500 bp in length
    print enzyme_entry, len(targetdigseqs[enzyme_entry])     
    if type(enzyme_entry) is tuple: #if the entry is a tuple (double digest)
        
        if enzyme_entry[0] in targetdigseqs: #if the first enzyme is in targetdigseqs (if a contig had only a single enzyme's cut sites)
            targetdigseqs[enzyme_entry] = targetdigseqs[enzyme_entry] + targetdigseqs[enzyme_entry[0]] #then add the single digest fragments to the double digest fragments
        if enzyme_entry[1] in targetdigseqs: #if the second enzyme is in feature digest (if a contig had only a single enzyme's cut sites)
            targetdigseqs[enzyme_entry] = targetdigseqs[enzyme_entry] + targetdigseqs[enzyme_entry[1]]  #then add the single digest fragments to the double digest fragments
        doubletargetdigseqs[enzyme_entry] = targetdigseqs[enzyme_entry]
 
print "doubletargetdigseqs"
for e in doubletargetdigseqs:
    print e, len(doubletargetdigseqs[e])
    
##NOTE 10-1-14 PM - featuredigest is not getting anything for single enzymes even with contig up to 5k, trying it with all_sites sorted (which should be case for these anyway)? Try to track down the cases where there are single cut sites and ask if there are any cases where the contig is in the gff file
#        
###write each unique sequence between 200 and 500 bp in length to a fasta file
#
##for enzyme in doubletargetdigseqs:
##    
##    targetseqfile = str(enzyme[0]) + "_" + str(enzyme[1]) + "_target_len_200_500_digestion_fragments.fasta"
##    description_string = str(enzyme[0]) + "_" + str(enzyme[1]) +"_digest_lettuce_contigs"
##    seqstowrite = [SeqRecord(Seq(o, generic_dna), id=str(n), description=description_string) for n,o in enumerate(set(targetdigseqs[enzyme]))]
##    SeqIO.write(seqstowrite, targetseqfile, "fasta")
##    
##
##################################################################################
#####  Categorize useful fragments by their location relative to genome features
#####  and write each enzyme's fragment lengths and C content to file
##################################################################################
##
##
##print "Begining section to write each enzyme's fragment lengths to file for generating a graphic."
##
##                   
##with open('RE_all_fragment_raw_data.csv', "wb") as csvfile:
##    summary = csv.writer(csvfile)
##    summary.writerow(["Enzyme", "Feature", "Proximity", "Length", "NumC"])
##         
##    for gff in doublefeaturedigest: # for gff file in feature digest
##        featurename = gff[:gff.find(".gff")]
##        for enz in doublefeaturedigest[gff]: #for enzyme (a dict key) in doublefeaturedigest[gff file]
##            print "Enzyme pair", enz
##            cat_counter = 0
##            for catlist in doublefeaturedigest[gff][enz]: #for a each list of lists containing [count C in frag, len frag, [5' cut site, 3' cut site]] for each fragment in a feature, the next for [count C in frag, len frag, [5' cut site, 3' cut site]] for fragment not in feature but within 1kb, the next for [count C in frag, len frag, [5' cut site, 3' cut site]] for fragment not within 1kb but within 5kb, the next not within 5kb but within 10kb, the next not within 10kb
##                #print len(catlist)                
##                if cat_counter ==0: #first time through loop (tracked by cut_counter) we're in the first list of doublefeaturedigest[gff][e], the list for in feature cuts
##                    prox = "within_feature"
##                elif cat_counter ==1: #second time through loop (tracked by cut_counter) we're in the second list of doublefeaturedigest[gff][e], the list for within 1kb of feature cuts
##                    prox = "proximal_1kb"
##                elif cat_counter ==2:  #third time through loop (tracked by cut_counter) we're in the third list of doublefeaturedigest[gff][e], the list for within 5kb of feature cuts
##                    prox = "proximal_5kb"
##                elif cat_counter ==3:  #fourth time through loop (tracked by cut_counter) we're in the fourth list of doublefeaturedigest[gff][e], the list for within 10kb of feature cuts
##                    prox = "proximal_10kb"
##                elif cat_counter ==4: #second time through loop (tracked by cut_counter) we're in the second list of doublefeaturedigest[gff][e], the list for beyond 10kb of feature cuts
##                    prox = "farther_than_10kb"  
##                
##                cat_counter +=1
##                
##                if catlist == []: #if for some reason one of the feature proximity bins is empty, skip it
##                    print "catlist for", prox, "is empty"
##                    continue
##                
##                for cut in catlist:
##                    #print cut
##                    summary.writerow([enz, gff, prox, cut[1], cut[0]])
##
##
