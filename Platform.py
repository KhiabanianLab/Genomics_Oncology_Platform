#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 11 11:57:59 2021

@author: nj277
"""
import tkinter as tk
from tkinter import filedialog, ttk, Canvas
from tkinter.filedialog import asksaveasfilename 
import pandas as pd
import numpy as np

# All-FIT and LOHGIC Libraries
from scipy import stats
from scipy.stats import binom
from scipy.stats import beta
import textwrap
import os
import xml.etree.ElementTree as ET

# --- classes ---
class MainWindow(tk.Tk): 
    def __init__(self):
        super().__init__()
        # root window
        self.geometry("1000x750")
        self.pack_propagate(False)
        self.resizable(0, 0)
        self.title('Genomics Oncology Platform')
        paddings = {'padx': 5, 'pady': 5}
        self.style = ttk.Style(self)
        self.style.theme_use('default')
        
        self.filename = None
        self.df = None
        
    
        # Title
        self.titleLabel = tk.Label(self,text = "Genomics Oncology Platform",font='Helvetica 18 bold')
        self.titleLabel.place(rely= 0.0, relx = 0.0)
        
        
        # Select File Type Drop Down Menu
        self.selectInputTypeDropDownLabel = tk.Label(self,text="Select input file type",**paddings,font='Helvetica 12 bold')
        self.selectInputTypeDropDownLabel.place(rely= 0.05, relx = 0.0)
        
        self.selectInputTypeDropDown_fileTypeVariable = tk.StringVar(self)
        self.selectInputTypeDropDown_fileTypeVariable.set("General File Format") # default value
        self.selectInputTypeDropDown = tk.OptionMenu(self, self.selectInputTypeDropDown_fileTypeVariable, "General File Format", "FoundationCDx xml File")
        self.selectInputTypeDropDown.place(rely= 0.05, relx = 0.13)
        self.selectInputTypeDropDown['width'] = 20
        
        # Load Data Button
        self.LoadFileButton = tk.Button(self, text='Load File',font='Helvetica 11 bold',command=self.load)
        self.LoadFileButton.place(rely= 0.1, relx = 0.01)  
        
        
        self.LoadFilePathLabel = tk.Label(self,bg="lightgray")
        self.LoadFilePathLabel.place(rely= 0.1, relx = 0.13,relwidth=0.83)
        
        
        # Treeview widget
        self.tv1 = ttk.Treeview(self)
        self.tv1.place(rely= 0.3, relx = 0.01,relheight=0.6,relwidth=0.95)
        self.treescrolly = tk.Scrollbar(self.tv1, orient="vertical", command=self.tv1.yview) # command means update the yaxis view of the widget
        self.treescrollx = tk.Scrollbar(self.tv1, orient="horizontal", command=self.tv1.xview) # command means update the xaxis view of the widget
        self.tv1.configure(xscrollcommand=self.treescrollx.set, yscrollcommand=self.treescrolly.set) # assign the scrollbars to the Treeview Widget
        self.treescrollx.pack(side="bottom", fill="x") # make the scrollbar fill the x axis of the Treeview widget
        self.treescrolly.pack(side="right", fill="y") # make the scrollbar fill the y axis of the Treeview widget
        
        # Biomarkers if present
        self.BiomarkersGroupLabel = tk.Label(self,text="Biomarkers and Purities",font='Helvetica 12 bold')
        self.BiomarkersGroupLabel.place(rely= 0.15, relx = 0.0)
        self.BiomarkersGroup =  tk.PanedWindow(self, relief = "groove", orient="horizontal")
        self.BiomarkersGroup.place(rely= 0.18, relx = 0.01,relheight=0.09,relwidth=0.95)
        
        self.TMB_score_status = tk.Label(self.BiomarkersGroup, text="TMB: ",font='Helvetica 11 bold')
        self.TMB_score_status.place(rely= 0.11, relx = 0.01)
        self.TMB_score_status_Value = tk.Label(self.BiomarkersGroup,bg="lightgray")
        self.TMB_score_status_Value.place(rely= 0.1, relx = 0.12,relwidth=0.1)
        
        self.MS_status = tk.Label(self.BiomarkersGroup, text="MSS: ",font='Helvetica 11 bold')
        self.MS_status.place(rely= 0.5, relx = 0.01)
        self.MS_status_Value = tk.Label(self.BiomarkersGroup,bg="lightgray")
        self.MS_status_Value.place(rely= 0.5, relx = 0.12,relwidth=0.1)
        
        # Purities
        self.PathologicalPurityLabel = tk.Label(self.BiomarkersGroup, text="Pathological Purity: ",font='Helvetica 11 bold')
        self.PathologicalPurityLabel.place(rely= 0.11, relx = 0.25)
        self.PathologicalPurity_Value = tk.Label(self.BiomarkersGroup,bg="lightgray")
        self.PathologicalPurity_Value.place(rely= 0.1, relx = 0.42,relwidth=0.1)
        
        self.ComputationalPurityLabel = tk.Label(self.BiomarkersGroup, text="Computational Purity: ",font='Helvetica 11 bold')
        self.ComputationalPurityLabel.place(rely= 0.5, relx = 0.25)
        self.ComputationalPurity_Value = tk.Label(self.BiomarkersGroup,bg="lightgray")
        self.ComputationalPurity_Value.place(rely= 0.5, relx = 0.42,relwidth=0.1)
        
        self.AllFITPurityLabel = tk.Label(self.BiomarkersGroup, text="All-FIT Purity: ",font='Helvetica 11 bold')
        self.AllFITPurityLabel.place(rely= 0.11, relx = 0.6)
        self.AllFITPurity_Value = tk.Label(self.BiomarkersGroup,bg="lightgray")
        self.AllFITPurity_Value.place(rely= 0.1, relx = 0.73,relwidth=0.1)
        
        
        # Clear Data Button
        self.ClearDataButton = tk.Button(self, text='Clear Data',font='Helvetica 11 bold', command=self.clearData)
        self.ClearDataButton.place(rely=0.91, relx= 0.01)
        
        # Save Data Button
        self.saveDataButton = tk.Button(self, text='Save Data',font='Helvetica 11 bold', command=self.saveData)
        self.saveDataButton.place(rely=0.91, relx= 0.86)
        
      
    def load(self):
        # check drop down selection
        if self.selectInputTypeDropDown_fileTypeVariable.get() == "General File Format":
            filename = filedialog.askopenfilename(filetypes=[('Python Files', '*.xlsx')])
        else:
            filename = filedialog.askopenfilename(filetypes=[('Python Files', '*.xml')])
        self.LoadFilePathLabel.config(text=filename)
        
        
        # Display Data 
        if filename:
            if filename.endswith('.xlsx'):
                self.df = pd.read_excel(filename)
            else:
                
                # get elements from Foundation xml and create a df
                tree = ET.parse(filename)
                tree_root = tree.getroot()
                
                Sample_ID = tree_root.find(".//TRFNumber").text
                Diagnosis = tree_root.find(".//SubmittedDiagnosis").text
                SpecimenSite = tree_root.find(".//SpecSite").text
                DOB = tree_root.find(".//DOB").text
                CollDate = tree_root.find(".//CollDate").text
                ReceivedDate = tree_root.find(".//ReceivedDate").text
                Gender = tree_root.find(".//Gender").text
                
                # Bait Set
                for sample in tree_root.iter('{http://foundationmedicine.com/compbio/variant-report-external}sample'):
                    baitSet = sample.get('bait-set')
                    meanExonDepth = sample.get('mean-exon-depth')
                
                # Tumor Mutation Burder
                TMB_all = TMBscore = TMBunit = TMBmutPerMb = TMBstatus = "unknown"
                for TMB in tree_root.iter('{http://foundationmedicine.com/compbio/variant-report-external}tumor-mutation-burden'):
                    if TMB is None:
                        continue
                    TMBscore = TMB.get('score')
                    TMBunit = TMB.get('unit')
                    TMBmutPerMb = TMB.get('mutations-per-megabase')
                    TMBstatus = TMB.get('status')
                    TMB_all = str(TMBscore)+"-"+TMBstatus
                self.TMB_score_status_Value.config(text=TMB_all,font='Helvetica 9 bold')
                    
                    
                    
                # Microsattelite instability status
                MSIstatus = "unknown"
                for MSI in tree_root.iter('{http://foundationmedicine.com/compbio/variant-report-external}microsatellite-instability'):
                    if MSI is None:
                        continue
                    MSIstatus = MSI.get('status')
                self.MS_status_Value.config(text=MSIstatus,font='Helvetica 9 bold')
                    
                
                # get pathological and computational purity
                for Var in tree_root.iter('{http://foundationmedicine.com/compbio/variant-report-external}variant-report'):
                    Computational_Purity = Var.get('purity-assessment')
                    Pathological_Purity  = Var.get('percent-tumor-nuclei')
                    
                    
                        
                    
                    
                
                # short variants
                Positions = []
                ShrVarGenes = []
                AAVars = []
                ShrVarMuts = []
                splitMutation = []
                ptMutations = []
                AllFreqs = []
                Depths = []
                Strands = []
                VariantClassification = []
                
                for shortVar in tree_root.iter('{http://foundationmedicine.com/compbio/variant-report-external}short-variant'):
                    position = shortVar.get('position')
                    Positions.append(position)
        
                    gene = shortVar.get('gene')
                    ShrVarGenes.append(gene)
        
                    aminoAcidVariant = shortVar.get('protein-effect')
                    AAVars.append(aminoAcidVariant)
        
                    alleleFreq = shortVar.get('allele-fraction')
                    alleleFreq = format(float(alleleFreq) * 100, '.2f')
                    AllFreqs.append(alleleFreq)
        
                    depth = shortVar.get('depth')
                    Depths.append(depth)
        
                    varClass = shortVar.get('functional-effect')
                    # rename variant classifications to proper maf identified classes
                    if varClass == "missense":
                        pvarClass = "Missense_Mutation"
                    elif varClass == "nonsense":
                        pvarClass = "Nonsense_Mutation"
                    elif varClass == "frameshift":
                    # check if ins or del
                        if "ins" in mutation:
                            pvarClass = "Frame_Shift_Ins"
                        elif "del" in mutation:
                            pvarClass = "Frame_Shift_Del"
                        else:
                            pvarClass = "Frame_Shift_Ins"
                    elif varClass == "nonframeshift":
                        if "ins" in mutation:
                            pvarClass = "In_Frame_Ins"
                        elif "del" in mutation:
                            pvarClass = "In_Frame_Del"
                        else:
                            pvarClass = "In_Frame_Ins"
                    else:
                        pvarClass = "Splice_Site"
        
                    VariantClassification.append(pvarClass)
        
                    strand = shortVar.get('strand')
                    Strands.append(strand)
                    mutation = shortVar.get('cds-effect')
                    if '>' in mutation:
                        splitmut = mutation.split(">")
                        ptMut = splitmut[0]
                        mutTo = splitmut[1]
                        truePtMut = ptMut[-1]
                        trueMutTo = mutTo
                    else:
                        splitmut = mutation
                        truePtMut = mutation
                        trueMutTo=""
            
                    splitMutation.append(splitmut)
                    ptMutations.append(truePtMut+" To "+trueMutTo)
                    ShrVarMuts.append(mutation)
                    
                # Copy number alterations
                CNA_genes = []
                CNA_copyNumbers = []
                CNA_types = []
                
                for CNA in tree_root.iter('{http://foundationmedicine.com/compbio/variant-report-external}copy-number-alteration'):
                    gene = CNA.get('gene')
                    CNA_genes.append(gene)
                    copyNumber = CNA.get('copy-number')
                    CNA_copyNumbers.append(copyNumber)
                    typeCNA = CNA.get('type')
                    CNA_types.append(typeCNA)
                    
                # Determine Copy Number (Ploidy) based on Overlap between Short Variants and Copy Number Alterations
                CNV = [2] * len(ShrVarGenes)
                intersection = [i for i in ShrVarGenes if i in CNA_genes];
                if not intersection:
                    print('no intersections')
                else:
                    for i in range(0,len(intersection)):
                        indx_shrVar = ShrVarGenes.index(''.join(intersection[i]))
                        indx_CNA = CNA_genes.index(''.join(intersection[i]))
                        cn = CNA_copyNumbers[indx_CNA]
                        CNV[indx_shrVar] = cn
                
                
                if Pathological_Purity is None:
                    self.df = pd.DataFrame(list(zip([Sample_ID] * len(ShrVarGenes),[Diagnosis] * len(ShrVarGenes),
                                                [SpecimenSite] * len(ShrVarGenes),[DOB] * len(ShrVarGenes),
                                                [CollDate] * len(ShrVarGenes),[ReceivedDate] * len(ShrVarGenes),
                                                [Gender] * len(ShrVarGenes), [TMB_all]*len(ShrVarGenes), [MSIstatus]*len(ShrVarGenes),
                                                ShrVarGenes, AAVars,Positions,ShrVarMuts,AllFreqs,
                                                Depths,CNV,VariantClassification,
                                                [Computational_Purity] * len(ShrVarGenes))),columns =['Sample_ID','Diagnosis','SpecimenSite','DOB',
                                                                                'CollectionDate','RecievedDate', 'Gender','TMB','MSI',
                                                                                'Gene', 'Protein_Change','Position','Mutation','VAF',
                                                                                'Depth','Copy_Number','Variant_Classification',
                                                                                'Computational_Purity'])
                elif Computational_Purity is None:
                    self.df = pd.DataFrame(list(zip([Sample_ID] * len(ShrVarGenes),[Diagnosis] * len(ShrVarGenes),
                                                [SpecimenSite] * len(ShrVarGenes),[DOB] * len(ShrVarGenes),
                                                [CollDate] * len(ShrVarGenes),[ReceivedDate] * len(ShrVarGenes),
                                                [Gender] * len(ShrVarGenes), [TMB_all]*len(ShrVarGenes), [MSIstatus]*len(ShrVarGenes),
                                                ShrVarGenes, AAVars,Positions,ShrVarMuts,AllFreqs,
                                                Depths,CNV,VariantClassification,[Pathological_Purity] * len(ShrVarGenes))),columns =['Sample_ID','Diagnosis','SpecimenSite','DOB',
                                                                                'CollectionDate','RecievedDate', 'Gender','TMB','MSI',
                                                                                'Gene', 'Protein_Change','Position','Mutation','VAF',
                                                                                'Depth','Copy_Number','Variant_Classification','Pathological_Purity'])
                else:                                                                                                                      
                # build dataframe
                    self.df = pd.DataFrame(list(zip([Sample_ID] * len(ShrVarGenes),[Diagnosis] * len(ShrVarGenes),
                                                [SpecimenSite] * len(ShrVarGenes),[DOB] * len(ShrVarGenes),
                                                [CollDate] * len(ShrVarGenes),[ReceivedDate] * len(ShrVarGenes),
                                                [Gender] * len(ShrVarGenes), [TMB_all]*len(ShrVarGenes), [MSIstatus]*len(ShrVarGenes),
                                                ShrVarGenes, AAVars,Positions,ShrVarMuts,AllFreqs,
                                                Depths,CNV,VariantClassification,[Pathological_Purity] * len(ShrVarGenes),
                                                [Computational_Purity] * len(ShrVarGenes))),columns =['Sample_ID','Diagnosis','SpecimenSite','DOB',
                                                                                'CollectionDate','RecievedDate', 'Gender','TMB','MSI',
                                                                                'Gene', 'Protein_Change','Position','Mutation','VAF',
                                                                                'Depth','Copy_Number','Variant_Classification','Pathological_Purity',
                                                                                'Computational_Purity'])
                
            
            self.filename = filename
            
            
            
            # Run All-FIT and LOHGIC
            All_FIT_Purity = self.RunAllFIT_LOHGIC()
            
            
            # add All_FIT_Purity to df
            self.df['All_FIT_Purity'] = All_FIT_Purity
            af = str((self.df['All_FIT_Purity'].to_list())[0])
            self.AllFITPurity_Value.config(text=af,font='Helvetica 9 bold')
            # get predicted model by purity
            
            ## check if both pathological and computational purities are present
            selectedModelCI = float(0)
            IDs = self.df['Sample_ID'].to_list()
            Genes = self.df['Gene'].to_list()
            VAFs = self.df['VAF'].to_list()
            Depths = self.df['Depth'].to_list()
            CNs = self.df['Copy_Number'].to_list()
            
            if 'Pathological_Purity' in self.df.columns:
                print("pathological purity present")
                pp = str((self.df['Pathological_Purity'].to_list())[0])
                self.PathologicalPurity_Value.config(text=pp,font='Helvetica 9 bold')
                # get models for pathological purity
                selectedPurity = self.df['Pathological_Purity'].to_list()
                bestModels,bestModel_weight,zygosity_list,loh_list = predictModelByPurity(selectedPurity,selectedModelCI,IDs,Genes,VAFs,Depths,CNs)
                # printed output: Best Model + W; Zygosity; LOH
                self.df['Pathological_Purity_Best_Model'] = [i +";"+"w="+ str(j) for i, j in zip(bestModels, bestModel_weight)]
                self.df['Pathological_Purity_Zygosity'] = zygosity_list
                self.df['Pathological_Purity_LOH'] = loh_list
            if 'Computational_Purity' in self.df.columns:
                print('computational purity present')
                cp = str((self.df['Computational_Purity'].to_list())[0])
                self.ComputationalPurity_Value.config(text=cp,font='Helvetica 9 bold')
                selectedPurity = self.df['Computational_Purity'].to_list()
                
                # get models for computational purity
                bestModels,bestModel_weight,zygosity_list,loh_list = predictModelByPurity(selectedPurity,selectedModelCI,IDs,Genes,VAFs,Depths,CNs)
                self.df['Computational_Purity_Best_Model'] = [i +";"+"w="+ str(j) for i, j in zip(bestModels, bestModel_weight)]
                self.df['Computational_Purity_Zygosity'] = zygosity_list
                self.df['Computational_Purity_LOH'] = loh_list
            
            ## get models for All-FIT purity
            selectedPurity = self.df['All_FIT_Purity'].to_list()
            bestModels,bestModel_weight,zygosity_list,loh_list = predictModelByPurity(selectedPurity,selectedModelCI,IDs,Genes,VAFs,Depths,CNs)
            self.df['All_FIT_Purity_Best_Model'] = [i +";"+"w="+ str(j) for i, j in zip(bestModels, bestModel_weight)]
            self.df['All_FIT_Purity_Zygosity'] = zygosity_list
            self.df['All_FIT_Purity_LOH'] = loh_list
                
            
            # display directly
            
            self.tv1["column"] = list(self.df.columns)
            self.tv1["show"] = "headings"
            for column in self.tv1["columns"]:
                self.tv1.heading(column, text=column) # let the column heading = column name
                
            self.df_rows = self.df.to_numpy().tolist() # turns the dataframe into a list of lists
            for row in self.df_rows:
                self.tv1.insert("", "end", values=row) # inserts each list into the treeview.
                
    def clearData(self):
        self.tv1.delete(*self.tv1.get_children())
        
    def saveData(self):
        
        savefile = filedialog.asksaveasfilename(defaultextension='.xlsx')
        
        self.df.to_excel(savefile, index=False, sheet_name="Results")  
        self.message = "Complete"
        #self.df.to_excel("output.xlsx")
        #self.label_text.set(self.message)
        
    def RunAllFIT_LOHGIC(self):
        
        
        # get unique sample IDs
        unique_Sample_ID = self.df.Sample_ID.unique()
        
        All_FIT_Purity =[]
        for sample in unique_Sample_ID:
            # get all rows for current sample
            temp_df = self.df[self.df['Sample_ID'] == sample]
            
            # set parameters
            std_dev = 2
            LOH_thres = 0.5
            W_thres = 0.7
            varType = 'all'
            CCF_pred_purity = ""
            CI_pred_purity = ""
            
            ind_pathologicalPurity = []
            ind_computationalPurity = []
            out_name = "tabledata_out"
            out_dir = os.getcwd()
            
            # get gene, vaf, depth and ploidy
            ind_SNV = temp_df['Gene'].tolist()
            ind_vaf = temp_df['VAF'].tolist()
            ind_depth = temp_df['Depth'].to_numpy()
            ind_ploidy = temp_df['Copy_Number'].to_numpy()
            
            
            CCF_pred_purity,CI_pred_purity = predicted_purity_from_CFF(ind_vaf,ind_depth,ind_SNV,out_dir+"/"+out_name,ind_ploidy,LOH_thres,std_dev,W_thres,varType)
            
            temp_AllFIT_Purity = [round(float(CCF_pred_purity)*100,2)] * len(ind_SNV) 
            All_FIT_Purity.extend(temp_AllFIT_Purity)
            
        return All_FIT_Purity
            
def predicted_purity_from_CFF(vaf_list,depth_list,SNV_list,out_name,ploidy_list,LOH_thres,std_dev,W_thres,varType):
    
    vaf_list = [float(i)/100 for i in vaf_list]
    depth_list = [float(i) for i in depth_list]
    
    all_purity = np.arange(0.01,1,0.01)
    CCF = [[] for x in range(len(all_purity))]
    weight = [[] for x in range(len(all_purity))]
    out = open(out_name+".txt","w")
    
    for pur in range(len(all_purity)):
        CCF[pur],weight[pur] = gray_box_accuracy(vaf_list,depth_list,all_purity[pur],ploidy_list,varType)
    
    gene_CCF_weight = {}
    for var in range(len(SNV_list)):
        gene_CCF_weight[SNV_list[var]] = []
    
    sum_CCF_weight = [0.0 for x in range(len(all_purity))]	#after removing germline no LOH and subclonal
    sum_CCF_weight_sq = [0.0 for x in range(len(all_purity))]
    Bef1_sum_CCF_weight = [0.0 for x in range(len(all_purity))]	#before removing any mutation
    Bef1_sum_CCF_weight_sq = [0.0 for x in range(len(all_purity))]
    Bef2_sum_CCF_weight = [0.0 for x in range(len(all_purity))]	#after removing germline no LOH
    Bef2_sum_CCF_weight_sq = [0.0 for x in range(len(all_purity))]
    
    for var in range(len(SNV_list)):
        name_likelihood = weight[0][var].keys()
        
        for pur in range(len(all_purity)):
            temp = 0.0
            for name in name_likelihood:
                temp += ((CCF[pur][var][name]-1)**2)*weight[pur][var][name]
            gene_CCF_weight[SNV_list[var]].append(np.log10(temp))
            Bef1_sum_CCF_weight[pur]+=temp
            Bef1_sum_CCF_weight_sq[pur] += temp**2
        
    ##remove germline mutation no LOH
    index_purity_bef1 = np.argmin(Bef1_sum_CCF_weight)
    if varType == "somatic":
        n_SNV, n_vaf, n_depth, n_ploidy = SNV_list, vaf_list, depth_list, ploidy_list
        Bef2_sum_CCF_weight = Bef1_sum_CCF_weight
        Bef2_sum_CCF_weight_sq = Bef1_sum_CCF_weight_sq
    else:
        germ_mut_index = []
        for mut in range(len(SNV_list)):
            germ_weight = 0.0
            flag_LOH = False
            if ploidy_list[mut] == 2:
                germ_list = ["Germline, LOH CNmut=1"]
            else:
                germ_list = []
            for i in range(ploidy_list[mut]):
                germ_list.append("Germline, CNmut=%i"%(i+1))
            
            for model in germ_list:
                germ_weight += weight[index_purity_bef1][mut][model]
                if model == "Germline, LOH CNmut=1" or model == "Germline, CNmut=%i"%ploidy_list[mut]:
                    if weight[index_purity_bef1][mut][model] > LOH_thres:
                        flag_LOH = True
            if not flag_LOH and germ_weight > W_thres:
                germ_mut_index.append(mut)
            
        
        n_SNV, n_vaf, n_depth, n_ploidy = remove_mut(vaf_list,depth_list,SNV_list,ploidy_list,germ_mut_index)
        for new_var in n_SNV:
            unlog = np.power(10,gene_CCF_weight[new_var])
            for pur in range(len(all_purity)):
                Bef2_sum_CCF_weight[pur] += unlog[pur]
                Bef2_sum_CCF_weight_sq[pur] += unlog[pur]**2
                
    ###remove subclonal mutation
    index_purity_bef2 = np.argmin(Bef2_sum_CCF_weight)
    p = all_purity[index_purity_bef2]
    subclonal_mut_index = []
    for mut in range(len(n_SNV)):
        f = n_vaf[mut]
        d = n_depth[mut]
        Y = n_ploidy[mut]
        if binom.cdf(round(d*f),d,p/(2*(1-p)+Y*p)) < 0.01:      #cumulative probability
            subclonal_mut_index.append(mut)
        
        n2_SNV, n2_vaf, n2_depth, n2_ploidy = remove_mut(n_vaf,n_depth,n_SNV,n_ploidy,subclonal_mut_index)
    
    for new_var in n2_SNV:
        unlog = np.power(10,gene_CCF_weight[new_var])
        for pur in range(len(all_purity)):
            sum_CCF_weight[pur] += unlog[pur]
            sum_CCF_weight_sq[pur] += unlog[pur]**2
            
    Bef1_sum_CCF_weight_LB = conf_interval(Bef1_sum_CCF_weight,Bef1_sum_CCF_weight_sq,len(SNV_list),std_dev)
    CI_p1 = all_purity[Bef1_sum_CCF_weight_LB < np.amin(Bef1_sum_CCF_weight)].tolist()
    CI_p2 = all_purity[conf_interval(Bef2_sum_CCF_weight,Bef2_sum_CCF_weight_sq,len(n_SNV),std_dev) < np.amin(Bef2_sum_CCF_weight)].tolist()
    sum_CCF_weight_LB = conf_interval(sum_CCF_weight,sum_CCF_weight_sq,len(n2_SNV),std_dev)
    CI_p3 = all_purity[sum_CCF_weight_LB < np.amin(sum_CCF_weight)].tolist()
    
    #print("purity_after_removing_germline_no_LOH_and_subclonal:\t"+str(round(all_purity[np.argmin(sum_CCF_weight)],3)))
            
        
    
    
    
    CCF_pred_purity = all_purity[np.argmin(sum_CCF_weight)]
    CI_pred_purity = 1
    return CCF_pred_purity,CI_pred_purity
       
def gray_box_accuracy(vaf_list,depth_list,purity,ploidy_list,varType):
    vaf_CI = 0.1
    f_incre = 0.005
    # vaf_list = [float(i)/100 for i in vaf_list]
    # depth_list = [float(i) for i in depth_list]
    
    CCF = [{} for j in range(len(vaf_list))]
    weight = [{} for j in range(len(vaf_list))]
    alpha = 1 - vaf_CI
    for j in range(len(vaf_list)):
        d = depth_list[j] 
        ploidy = ploidy_list[j]
        if ploidy == 2:
            if varType == "somatic":
                name_likelihood = ["Somatic, LOH CNmut=1"]
            else:
                name_likelihood = ["Somatic, LOH CNmut=1","Germline, LOH CNmut=1"]
        else:
            name_likelihood = []
        
        som_list = []
        germ_list = []
        for i in range(ploidy):
            som_list.append("Somatic, CNmut=%i"%(i+1))
            germ_list.append("Germline, CNmut=%i"%(i+1))
            
        name_likelihood[1:1] = som_list
        if varType != "somatic":
            name_likelihood.extend(germ_list)
        
        for name in name_likelihood:
            weight[j][name] = 0.0 # initialize weights
        
        if d == 0:
            continue
        if ploidy == 2:
            CCF[j]["Somatic, LOH CNmut=1"] = vaf_list[j]/(purity/(2-purity)) #somatic LOH
            if varType != "somatic":
                CCF[j]["Germline, LOH CNmut=1"] = ((vaf_list[j]*(2-purity))-1+purity)/purity	#germline LOH
        
        for i in range(ploidy):
            CCF[j]["Somatic, CNmut=%i"%(i+1)] = vaf_list[j]/(((i+1)*purity)/(2*(1-purity)+ploidy*purity))	#somatic LOH high CN
            if varType != "somatic":
                CCF[j]["Germline, CNmut=%i"%(i+1)] = ((vaf_list[j]*(2-2*purity+ploidy*purity))-1+purity)/((i+1)*purity)	#germline LOH high CN
        
        freqLB,freqUB = binom_interval(round(d*vaf_list[j]),d,alpha)
        
        f_range = np.arange(freqLB,freqUB+f_incre,f_incre)
        f_range[-1] = freqUB
        
        # compute aic weights
        aic = [{} for f in range(len(f_range))]
        for f in range(len(f_range)):
            freq = f_range[f]
            p = purity
            if ploidy == 2:		#dbinom-exact probability
                t_ans1 = binom.pmf(round(d*freq),d,p/(2-p)) + np.finfo(np.double).tiny     #somatic LOH
                aic[f]["Somatic, LOH CNmut=1"] = 2-2*np.log(t_ans1)
                if varType != "somatic":
                    t_ans2 = binom.pmf(round(d*freq),d,1/(2-p)) + np.finfo(np.double).tiny     #germline LOH
                    aic[f]["Germline, LOH CNmut=1"] = 2-2*np.log(t_ans2)
            
            for i in range(ploidy):
                t_ans3 = binom.pmf(round(d*freq),d,((i+1)*p)/(2*(1-p)+ploidy*p)) + np.finfo(np.double).tiny        #somatic LOH high CN
                aic[f]["Somatic, CNmut=%i"%(i+1)] = 2-2*np.log(t_ans3)
                if varType != "somatic":
                    t_ans4 = binom.pmf(round(d*freq),d,(1-p+(i+1)*p)/(2*(1-p)+ploidy*p)) + np.finfo(np.double).tiny    #germline LOH high CN
                    aic[f]["Germline, CNmut=%i"%(i+1)] = 2-2*np.log(t_ans4)
                
        smallest = aic[0]["Somatic, CNmut=1"]
        for f in range(len(f_range)):
            for name in name_likelihood:
                if smallest > aic[f][name]:
                    smallest = aic[f][name]
        
        D = 0.0
        for f in range(len(f_range)):
            for name in name_likelihood:
                D += np.exp(-0.5*(aic[f][name]-smallest))
                
        for f in range(len(f_range)):
            for name in name_likelihood:
                weight[j][name] += (np.exp(-0.5*(aic[f][name]-smallest))/D)
                
    return CCF, weight
    
def binom_interval(n_success, total, conf_int):
    
    quantile = (1 - conf_int)/2
    lower = beta.ppf(quantile, n_success, total - n_success + 1)
    upper = beta.ppf(1 - quantile, n_success + 1, total - n_success)
    if n_success == total:
        upper = 1.0
    if n_success == 0:
        lower = 0.0
    
    
    return (lower, upper)

def remove_mut(vaf_list,depth_list,SNV_list,ploidy_list,index):
	n_SNV = []
	n_vaf = []
	n_depth = []
	n_ploidy = []
	for ind in range(len(SNV_list)):
		if ind not in index:
			n_SNV.append(SNV_list[ind])
			n_vaf.append(vaf_list[ind])
			n_depth.append(depth_list[ind])
			n_ploidy.append(ploidy_list[ind])
	return n_SNV, n_vaf, n_depth, n_ploidy

def conf_interval(x,x2,num_mut,std_dev):
	x_err = ((np.array(x2)/num_mut) - (np.array(x)/num_mut)**2)**0.5 
	return x-std_dev*np.array(x_err)


def predictModelByPurity(selectedPurity,selectedModelCI,IDs,Genes,VAFs,Depths,CNs):
    
    selectedPurity = [float(i)/100 for i in selectedPurity]
    
    
    VAFs = [float(i)/100 for i in VAFs]
    Depths = [int(i) for i in Depths]
    CNs = [int(i) for i in CNs]
    # set parameters 
    vaf_CI = 0.1
    f_incre = 0.005
    p_incre = 0.01
    alpha = 1 - vaf_CI
    LOH_thres = 0.5
    W_thres = 0.7
    
    # initialize
    bestModels = []
    bestModel_weight = []
    zygosity_list = []
    loh_list = []
    
    for j in range(len(Genes)):
        weight = {}
        d = Depths[j] 
        ploidy = CNs[j]
        purity = selectedPurity[j]
        
        if ploidy == 2:
            name_likelihood = ["Somatic LOH CNmut=1","Germline LOH CNmut=1"]
        else:
            name_likelihood = ["Somatic no LOH","Germline no LOH"] #added
        
        som_list = []
        germ_list = []
        
        for i in range(ploidy):
            som_list.append("Somatic CNmut=%i"%(i+1))
            germ_list.append("Germline CNmut=%i"%(i+1))
        
        name_likelihood[1:1] = som_list
        name_likelihood.extend(germ_list)
        
        for name in name_likelihood:
            weight[name] = 0.0
        
        if d == 0:
            continue
        
        freqLB,freqUB = new_binom_interval(round(d*VAFs[j]),d,alpha)
        f_range = np.arange(freqLB,freqUB+f_incre,f_incre)
        f_range[-1] = freqUB
        
        # determine confidence interval
        if selectedModelCI == 0:
            CI = [(purity-0.05),(purity+0.05)]
        else:
            CI = [float(selectedModelCI[0]),float(selectedModelCI[len(selectedModelCI)-1])]
        
        
        if CI[0] == CI[1]:
            pur_range = np.array([CI[0]])
        else:
            pur_range = np.arange(CI[0],CI[1]+p_incre,p_incre)
            pur_range[-1] = CI[1]
        
        aic = [[{} for f in range(len(f_range))] for pur in range(len(pur_range))]
        for pur in range(len(pur_range)):
            for f in range(len(f_range)):
                freq = f_range[f]
                p = pur_range[pur]
                
                if ploidy == 2:		#dbinom-exact probability
                    t_ans1 = binom.pmf(round(d*freq),d,p/(2-p)) + np.finfo(np.double).tiny     #somatic LOH
                    aic[pur][f]["Somatic LOH CNmut=1"] = 2-2*np.log(t_ans1)
                    t_ans2 = binom.pmf(round(d*freq),d,1/(2-p)) + np.finfo(np.double).tiny     #germline LOH
                    aic[pur][f]["Germline LOH CNmut=1"] = 2-2*np.log(t_ans2)
                
                for i in range(ploidy):
                    t_ans3 = binom.pmf(round(d*freq),d,((i+1)*p)/(2*(1-p)+ploidy*p)) + np.finfo(np.double).tiny        #somatic LOH high CN
                    aic[pur][f]["Somatic CNmut=%i"%(i+1)] = 2-2*np.log(t_ans3)
                    t_ans4 = binom.pmf(round(d*freq),d,(1-p+(i+1)*p)/(2*(1-p)+ploidy*p)) + np.finfo(np.double).tiny    #germline LOH high CN
                    aic[pur][f]["Germline CNmut=%i"%(i+1)] = 2-2*np.log(t_ans4)
                
        smallest = aic[0][0]["Somatic CNmut=1"]
        for pur in range(len(pur_range)):
            for f in range(len(f_range)):
                for name in name_likelihood:
                    if smallest > aic[pur][f][name]:
                        smallest = aic[pur][f][name]
        
        D = 0.0
        for pur in range(len(pur_range)):
            for f in range(len(f_range)):
                for name in name_likelihood:
                    D += np.exp(-0.5*(aic[pur][f][name]-smallest))
        
        for pur in range(len(pur_range)):
            for f in range(len(f_range)):
                for name in name_likelihood:
                    weight[name] += (np.exp(-0.5*(aic[pur][f][name]-smallest))/D)
        
        best_model = ""
        largest_w = 0
        zygosity = ""
        LOH_status = ""
        
        for name in name_likelihood:
            if weight[name] > largest_w:
                largest_w = weight[name]
                best_model = name
        
                
        
        if binom.cdf(round(d*VAFs[j]),d,purity/(2*(1-purity)+ploidy*purity)) < 0.01:
            best_model = best_model+", subclonal"
        
        bestModels.append(best_model)
        bestModel_weight.append(round(largest_w,2))
        
        
        all_germ = germ_list[:]
        all_som = som_list[:]
        germ_w = 0.0
        som_w = 0.0
        
        if ploidy == 2:
            all_germ.insert(0,"Germline LOH CNmut=1")
            all_som.insert(0,"Somatic LOH CNmut=1")
        
        for model in all_germ:
            germ_w += weight[model]
        for model in all_som:
            som_w += weight[model]
        
        flag_LOH = 0
        if germ_w > som_w and germ_w > W_thres:	#germline
            zygosity = "Germline"
            for model in all_germ:
                if model == "Germline LOH CNmut=1" or model == "Germline CNmut=%i"%ploidy:
                    if weight[model] > LOH_thres:
                        flag_LOH = 1
        elif som_w > germ_w and som_w > W_thres:
            zygosity = "Somatic"
            for model in all_som:
                if model == "Somatic LOH CNmut=1" or model == "Somatic CNmut=%i"%ploidy:
                    if weight[model] > LOH_thres:
                        flag_LOH = 1
        else:
            zygosity = "Ambiguous"
            flag_LOH = -1
        
        if flag_LOH == 1:
            LOH_status = "Yes"
        elif flag_LOH == 0:
            LOH_status = "No"
        else:
            LOH_status = "Ambiguous"
        
        zygosity_list.append(zygosity)
        loh_list.append(LOH_status)
        
    
    
    return bestModels,bestModel_weight,zygosity_list,loh_list

def new_binom_interval(n_success, total, conf_int):
    quantile = (1 - conf_int) / 2
    
    lower = beta.ppf(quantile, n_success, total - n_success + 1)
    upper = beta.ppf(1 - quantile, n_success + 1, total - n_success)
    
    if n_success == total:
        upper = 1.0
        
    if n_success == 0:
        lower = 0.0
        
    return (lower, upper)
        

# --- main ---

if __name__ == '__main__':
    #root = tk.Tk()
    # root.geometry("500x500")
    # root.pack_propagate(False)
    # # Setting Theme
    # style = ThemedStyle(root)
    # style.set_theme("equilux")

    # root.configure(bg='white')
    
    app = MainWindow()
    app.mainloop()        
        
        
