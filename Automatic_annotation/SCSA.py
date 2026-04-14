#########################################################################
# File Name: scRNA_anno.py
# > Author: CaoYinghao
# > Mail: caoyinghao@gmail.com 
#########################################################################
#! /usr/bin/python

import sys
import argparse
import gzip
import os

import numpy as np

def get_final_path(filepath):
    """Helper function to ensure filepath is in the final folder"""
    # Create final directory if it doesn't exist
    final_dir = "final"
    if not os.path.exists(final_dir):
        os.makedirs(final_dir)
    
    # Get the directory and filename
    dirname = os.path.dirname(filepath) if os.path.dirname(filepath) else ""
    filename = os.path.basename(filepath)
    
    # If there's a subdirectory structure, preserve it in final/
    if dirname:
        final_path = os.path.join(final_dir, dirname, filename)
        # Create subdirectory if needed
        subdir = os.path.join(final_dir, dirname)
        if not os.path.exists(subdir):
            os.makedirs(subdir)
    else:
        final_path = os.path.join(final_dir, filename)
    
    return final_path
from numpy import asfarray,arange,minimum,abs,mat,sum,power,mean,array,std,log2,log10
import pandas as pd
from pandas import DataFrame,read_csv,ExcelWriter
from pickle import dump,load
from scipy.stats import fisher_exact
from scipy.sparse import coo_matrix

class Annotator(object):
    def __init__(self,args):
        self.args = args
        pass

    @staticmethod
    def do_fisher_test(x,num1,num2):
        """do fisher-exact test for go annotation"""
        fy = x['gene_num']
        by = x['othergene_num']
        o,p = fisher_exact([[fy,num1 - fy],[by,num2 - by]],alternative="greater")
        return p
        pass

    @staticmethod
    def do_sig_tag(x):
        if x <= 0.001:
            return "***"
        elif x <= 0.01:
            return "**"
        elif x <= 0.05:
            return "*"
        else:
            return "-"

    @staticmethod
    def p_adjust_bh(p):
        """Benjamini-Hochberg p-value correction for multiple hypothesis testing."""
        p = asfarray(p)
        by_descend = p.argsort()[::-1]
        by_orig = by_descend.argsort()
        steps = float(len(p)) / arange(len(p), 0, -1)
        q = minimum(1, minimum.accumulate(steps * p[by_descend]))
        return q[by_orig]

    @staticmethod
    def to_output(h_values,wb,outtag,cname,title):
        if outtag.lower() == "ms-excel":
            h_values.to_excel(wb,sheet_name = "Cluster " + cname + " " + title,index=False)
        else:
            # Convert all string columns to ASCII
            for col in h_values.columns:
                if h_values[col].dtype == 'object':
                    h_values[col] = h_values[col].astype(str).apply(lambda x: x.encode('ascii', 'ignore').decode('ascii'))
            
            # Handle both file paths and file handles
            if isinstance(wb, str):
                with open(wb, 'w', encoding='ascii', errors='ignore') as f:
                    h_values.to_csv(f, sep="\t", quotechar="\t", index=False, header=False)
            else:
                h_values.to_csv(wb, sep="\t", quotechar="\t", index=False, header=False)
        pass

    @staticmethod
    def translate_go(name="go.obo"):
        """extract go annotation dicts from obo"""
        ids = []
        names = []
        cls = []
        infile = None
        if name.endswith(".gz"):
            infile = gzip.open(name,"rt")
        else:
            infile = open(name,"rt")
        for line in infile:
            if line.strip().startswith("id:"):
                ids.append(line.strip()[4:])
            if line.strip().startswith("name:"):
                names.append(line.strip()[6:])
            if line.strip().startswith("namespace:"):
                cls.append(line.strip()[11:])
        infile.close()
        return dict(zip(ids,names))


    def do_go_annotation(self,gof,fore,back,cname,gtype):
        """return go annotation with significance tag"""
        fil = gof[2].map(lambda value: len(set([value]) & fore) > 0)
        fgnames = gof[fil].groupby(by=4)[2].unique()
        bfil = gof[2].map(lambda value: len(set([value]) & back) > 0)
        bgnames = gof[bfil].groupby(by=4)[2].unique()
        dat = DataFrame({"genes":fgnames,"othergenes":bgnames})
        num1 = len(fore)
        num2 = len(back)
        dat = dat.dropna(axis=0)
        dat['gene_num'] = dat['genes'].map(lambda x:x.size)
        dat['othergene_num'] = dat['othergenes'].map(lambda x:x.size)
        if dat.size ==0:
            return dat

        dat['p-value'] =  dat.apply(Annotator.do_fisher_test,axis=1,args=(num1,num2))
        used = dat.loc[:,dat.columns]   #.head(5)
        used['ids'] = used.index
        used['q-value'] = Annotator.p_adjust_bh(used['p-value'])
        used['sig'] = used['q-value'].map(Annotator.do_sig_tag)
        used['go_name'] = used['ids'].map(lambda x:self.gos[x])
        
        def convert_to_string(x):
            if isinstance(x, np.ndarray):
                # Handle numpy arrays (including object dtype arrays)
                return '; '.join(str(g) for g in x)
            elif isinstance(x, (list, tuple, set)):
                return '; '.join(str(g) for g in x)
            elif isinstance(x, str):
                return x
            else:
                return str(x)
        
        used['genes_in_cluster'] = used['genes'].map(convert_to_string).astype(str)
        used['genes_outside_cluster'] = used['othergenes'].map(convert_to_string).astype(str)
        used['cluster'] = cname
        outs = used.sort_values(by='p-value')[['ids','cluster','gene_num','othergene_num','p-value','q-value',"sig",'go_name','genes_in_cluster','genes_outside_cluster']]
        if self.args.output:
            # Replace single letter codes with full descriptive names
            go_class_mapping = {
                'P': 'Biological Process',
                'F': 'Molecular Function', 
                'C': 'Cellular Component'
            }
            outs['go_class'] = go_class_mapping.get(gtype, gtype)
        return outs
        pass

    def print_class(self,h_values,cname):
        """print cell predictions with scores."""
        o = ""
        titlebar = "-" * 60 + "\n"
        
        # Save top 5 cell types and scores to Excel

        if h_values is not None and hasattr(h_values, 'size') and h_values.size > 0:
            top5 = h_values.head(5)
            cluster_num = cname
            result_df = pd.DataFrame({
                'Cluster': [cluster_num] * len(top5),
                'Cell Type': top5['Cell Type'].values,
                'Score': top5['Z-score'].values if 'Z-score' in top5.columns else top5.iloc[:, 1].values
            })
            input_basename = os.path.splitext(os.path.basename(self.args.input))[0]
            excel_file = f'{input_basename}_top5_celltype_prediction.xlsx'
            excel_file = get_final_path(excel_file)
            if os.path.exists(excel_file):
                existing_df = pd.read_excel(excel_file)
                result_df = pd.concat([existing_df, result_df], ignore_index=True)
            result_df.to_excel(excel_file, index=False)
            
            # Save Good or ? results (top 1 for Good, top 2 for ?)
            result_type = None
            good_or_question_df = None
            cluster_num = cname
            if h_values.size == 3 or h_values.size == 2:
                result_type = "Good"
                cell_type = h_values.values[0][0]
                score = h_values.values[0][1]
                good_or_question_df = pd.DataFrame({
                    'Cluster': [cluster_num],
                    'Cell Type': [cell_type],
                    'Score': [score]
                })
            elif h_values.size > 2:
                if float(h_values.iloc[0,1])/float(h_values.iloc[1,1]) >= 2 or float(h_values.iloc[1,1]) < 0:
                    result_type = "Good"
                    good_or_question_df = pd.DataFrame({
                        'Cluster': [cluster_num],
                        'Cell Type': [h_values['Cell Type'].values[0]],
                        'Score': [h_values['Z-score'].values[0]]
                    })
                else:
                    result_type = "?"
                    good_or_question_df = pd.DataFrame({
                        'Cluster': [cluster_num] * 2,
                        'Cell Type': [h_values['Cell Type'].values[0], h_values['Cell Type'].values[1]],
                        'Score': [h_values['Z-score'].values[0], h_values['Z-score'].values[1]]
                    })
            
            if good_or_question_df is not None:
                good_excel_file = f'{input_basename}_celltype_predictions.xlsx'
                good_excel_file = get_final_path(good_excel_file)
                if os.path.exists(good_excel_file):
                    existing_good_df = pd.read_excel(good_excel_file)
                    good_or_question_df = pd.concat([existing_good_df, good_or_question_df], ignore_index=True)
                good_or_question_df.to_excel(good_excel_file, index=False)
            
       
        if h_values is None:
            if self.args.noprint:
                return "E",None,"-","-","-"
            o += titlebar
            o += "{0:<10}{1:^30}{2:<10}".format("Type","Cell Type","Score")
            o += "\n" + "-"*60 + "\n"
            o += "{0:<10}{1:^30}{2:<10}".format("-","-","-")
            o += "\n" + titlebar
            return "E",None,"-","-","-"
        elif h_values.size == 0:
            if self.args.noprint:
                return "N",None,"-","-","-"
            o += titlebar
            o += "{0:<10}{1:^30}{2:<10}".format("Type","Cell Type","Score")
            o += "\n" + "-"*60 + "\n"
            o += "{0:<10}{1:^30}{2:<10}".format("-","-","-")
            o += "\n" + titlebar
            return "N",None,"-","-","-"
        elif h_values.size == 3:
            if self.args.noprint:
                return "Good",o,h_values.values[0][0],h_values.values[0][1],"-"
            o += titlebar
            o += "{0:<10}{1:^30}{2:<10}{3:<5}".format("Type","Cell Type","Score","Times")
            o += "\n" + "-"*60 + "\n"
            o += "{0:<10}{1:^30}{2:<10.4f}".format("Good",h_values.values[0][0],h_values.values[0][1])
            o += "\n" + titlebar
            return "Good",o,h_values.values[0][0],h_values.values[0][1],"-"
            pass
        elif h_values.size == 2:
            if self.args.noprint:
                return "Good",o,h_values.values[0][0],h_values.values[0][1],"-"
            o += titlebar
            o += "{0:<10}{1:^30}{2:<10}{3:<5}".format("Type","Cell Type","Score","Times")
            o += "\n" + "-"*60 + "\n"
            o += "{0:<10}{1:^30}{2:<10.4f}".format("Good",h_values.values[0][0],h_values.values[0][1])
            o += "\n" + titlebar
            return "Good",o,h_values.values[0][0],h_values.values[0][1],"-"
            pass
        elif float(h_values.iloc[0,1])/float(h_values.iloc[1,1]) >= 2 or float(h_values.iloc[1,1] < 0):
            times = np.abs(float(h_values.iloc[0,1])/float(h_values.iloc[1,1]))
            if self.args.noprint:
                return "Good",o,h_values.values[0][0],h_values.values[0][1],times
            o += titlebar
            o += "{0:<10}{1:^30}{2:<10}{3:<5}".format("Type","Cell Type","Score","Times")
            o += "\n" + titlebar
            o += "{0:<10}{1:^30}{2:<10.4f}{3:<5.1f}".format("Good",h_values['Cell Type'].values[0],h_values['Z-score'].values[0],times)
            o += "\n" + titlebar
            return "Good",o,h_values['Cell Type'].values[0],h_values['Z-score'].values[0],times
            pass
        else:
            times = np.abs(float(h_values.iloc[0,1])/float(h_values.iloc[1,1]))
            if self.args.noprint:
                return "?",o,str(h_values['Cell Type'].values[0]) + "|" + str(h_values['Cell Type'].values[1]),str(h_values['Z-score'].values[0]) + "|" + str(h_values['Z-score'].values[1]),times
            o += titlebar
            o += "{0:<10}{1:^30}{2:<10}{3:<5}".format("Type","Cell Type","Score","Times")
            o += "\n" + titlebar
            o += "{0:<10}{1:^30}{2:<10.4f}{3:<5.1f}".format("?",h_values['Cell Type'].values[0],h_values['Z-score'].values[0],times)
            o += "\n" + titlebar
            o += "{0:<10}{1:^29}({2:<.4f})".format("","("+h_values['Cell Type'].values[1]+")",h_values['Z-score'].values[1])
            o += "\n" + titlebar
            return "?",o,str(h_values['Cell Type'].values[0]) + "|" + str(h_values['Cell Type'].values[1]),str(h_values['Z-score'].values[0]) + "|" + str(h_values['Z-score'].values[1]),times
        pass

    def deal_with_badtype(self,cname,other_gene_names,colnames):
        """go annotation need to be performed"""
        if len(self.human_gofs) != 0:
            fset = set()
            bset = set()
            if colnames is None:
                print("!WARNING(go processing):Zero gene sets found for the cluster",cname)
                print("!WARNING(go processing):Change the threshold and try again?")
                return 
            for c in colnames:
                if c in self.ensem_hgncs:
                    fset.add(self.ensem_hgncs[c])
                else:
                    fset.add(c)
            for c in other_gene_names:
                if c in self.ensem_hgncs:
                    bset.add(self.ensem_hgncs[c])
                else:
                    bset.add(c)
            if len(fset) == 0:
                print("!WARNING(go processing):Zero gene sets found for the cluster",cname)
                print("!WARNING(go processing):Change the threshold and try again?")
                return 
            if len(bset) == 0:
                print("!WARNING(go processing):Zero gene sets found for other clusters")
                print("!WARNING(go processing):Change the threshold and try again?")
                return 
            names = ["Function","Component","Process"]
            if self.args.noprint == False:
                print("Go Enrichment analysis:","Group1:",len(fset),"Group2:",len(bset))
            if len(fset) > 0 and len(bset) > 0:
                all_outs = DataFrame()
                for i,f in enumerate(self.human_gofs):
                    o = " ".join([">"*30,names[i], "<"*30])
                    if self.args.noprint == False:
                        print(o)
                    outs = self.do_go_annotation(f,fset,bset,cname,names[i][0])
                    if outs.size == 0:continue
                    if all_outs.size == 0:
                        all_outs = outs
                    else:
                        all_outs = pd.concat([all_outs, outs], ignore_index=True)
                    if self.args.noprint == False:
                        print()
                if self.args.output:
                    # Store results for later concatenation instead of writing immediately
                    if not hasattr(self, 'all_go_results'):
                        self.all_go_results = []
                    self.all_go_results.append(all_outs)
                    print(f"Added GO results for cluster {cname}, total results: {len(self.all_go_results)}")



    def calcu_cellranger_group(self,expfile,hgvc=False):
        """deal with cellranger input matrix"""
        exps = read_csv(expfile)
        columns = exps.columns

        pre,suf,suf1 ="Cluster "," UMI counts/cell"," Weight"
        fid = "Gene Name" if hgvc == True else "Gene ID"
        gcol = "gene" if hgvc == True else "ensemblID"
        ccol = "cellName"

        if self.args.target.lower() not in ["cancersea","cellmarker"]:
            print("Error target : -t, --target,(cellmarker,[cancersea])")
            sys.exit(0)

        if self.args.target.lower() == "cancersea":
            gcol = "gene" if hgvc == True else "ensemblID"
            ccol = "name"


        abs_tag = False

        cnum = int(len(exps.columns) / 2 - 1)
        ver_tag = "V1"
        pname = ""

        if "Feature ID" in columns: # v3
            fid = "Feature Name" if hgvc == True else "Feature ID"
            pre,suf,suf1 ="Cluster "," Mean Counts"," Log2 fold change"
            cnum = int((len(exps.columns)-2) / 3)
            pname = " Adjusted p value"
            self.args.weight = self.args.foldchange
            ver_tag = "V3"
        elif "Cluster 1 Mean UMI Counts" in columns: # v2
            fid = "Gene Name" if hgvc == True else "Gene ID"
            pre,suf,suf1 ="Cluster "," Mean Counts"," Log2 fold change"
            cnum = int((len(exps.columns)-2) / 3)
            self.args.weight = self.args.foldchange
            pname = " Adjusted p value"
            ver_tag = "V2"
        
        # Extract cluster names from column headers
        cluster_names = set()
        for col in columns:
            if col.startswith(pre) and (suf1 in col or " Weight" in col):
                # Extract cluster name: "Cluster {name} {suffix}" -> "{name}"
                if ver_tag == "V1":
                    if " Weight" in col:
                        cluster_name = col.replace(pre, "").replace(" Weight", "")
                        cluster_names.add(cluster_name)
                elif ver_tag == "V2" or ver_tag == "V3":
                    if suf1 in col:
                        cluster_name = col.replace(pre, "").replace(suf1, "")
                        cluster_names.add(cluster_name)
        
        # If no named clusters found, fall back to numeric (backward compatibility)
        if not cluster_names:
            cluster_names = [str(i) for i in range(1, cnum+1)]
        else:
            # Sort cluster names: numeric first, then alphabetical
            def sort_key(name):
                try:
                    return (0, int(name))  # Numeric clusters first
                except ValueError:
                    return (1, name)  # Named clusters after
            cluster_names = sorted(cluster_names, key=sort_key)
        
        outs = []

        self.wb = self.wbgo = None
        if self.args.output:
            # Save to final folder
            original_output = self.args.output
            # self.args.output = get_final_path(self.args.output)
            if self.args.outfmt.lower() == "ms-excel":
                if not self.args.output.endswith(".xlsx") and (not self.args.output.endswith(".xls")):
                    self.args.output += ".xlsx"
                self.wb = ExcelWriter(self.args.output)
                self.wbgo = self.wb
            elif self.args.outfmt.lower() == "txt":
                self.wb = open(self.args.output,"w")
                self.wb.write("Cell Type\tZ-score\tCluster\n")
                self.wbgo = open(self.args.output + ".go","w")
                self.wbgo.write('ids\tgene_num\tothergene_num\tp-value\tq-value\tsig\tname\tcluster\tgo_class\n')
            else:
                print("Error output format: -m, --outfmt,(ms-excel,[txt])")
                sys.exit(0)

        

        for cname in cluster_names:
            if self.args.cluster != "all":
                if self.args.cluster.find(",") > -1:
                    sets = self.args.cluster.split(",")
                    if cname not in sets:
                        continue
                else:
                    if cname != self.args.cluster:
                        continue
            #if i != 1 :continue
            o = " ".join(["#"*30,"Cluster",cname, "#"*30]) + "\n"
            if self.args.noprint == False:
                print(o)
            ptitle = pre + cname + pname
            ltitle = pre + cname + suf1
            if ltitle not in exps.columns:
                print(ltitle,"column not in the input table!")
                sys.exit(0)

            newexps = None
            if ver_tag == "V1":
                newexps = exps[exps[ltitle]>=self.args.weight]
            else:
                newexps = exps[(exps[ltitle]>=self.args.weight) & (exps[ptitle] <= self.args.pvalue)]
            #print(newexps.shape)
            h_values,colnames = self.get_cell_matrix(newexps,ltitle,fid,gcol,ccol,abs_tag)
            if h_values is None:
                t,o_str,c,v,times = self.print_class(h_values,cname)
                outs.append([cname,t,c,v,times])
                if self.args.noprint == False:
                    print(o_str)
                continue
            h_values['Cluster'] = cname
            if self.args.output:
                Annotator.to_output(h_values,self.wb,self.args.outfmt,cname,"Cell Type")


            t,o_str,c,v,times = self.print_class(h_values,cname)
            outs.append([cname,t,c,v,times])
            if self.args.noprint == False:
                print(o_str)
            other_gene_names = set()
            for other_cname in cluster_names:
                if cname == other_cname: continue
                jtitle = pre + other_cname + suf1
                jptitle = pre + other_cname + pname
                otherexps = None
                if ver_tag == "V1":
                    otherexps = exps[exps[jtitle]>=self.args.weight]
                else:
                    otherexps = exps[(exps[jtitle]>=self.args.weight) & (exps[jptitle] <= self.args.pvalue)]
                if self.args.target.lower() == "cancersea":
                    tfc,trownames,trownum,tcolnames,tcolnum = self.get_cell_gene_names(otherexps,self.smarkers,fid,gcol,ccol,"other")
                    other_gene_names |= set(tcolnames)
                elif self.args.target.lower() == "cellmarker":
                    tfc,trownames,trownum,tcolnames,tcolnum = self.get_cell_gene_names(otherexps,self.cmarkers,fid,gcol,ccol,"other")
                    if not trownames:
                        #print("WARNING3:Zero gene sets found for the cluster" + str(j))
                        #print("WARNING3:Change the threshold and try again?")
                        continue
                    other_gene_names |= set(tcolnames)
            #print("Other Gene number:",len(other_gene_names))
            self.deal_with_badtype(cname,other_gene_names,colnames)
        if self.args.output:
            self.wb.close()
            self.wbgo.close()
            # Save combined GO results to single CSV file
            print(f"Checking for GO results: hasattr={hasattr(self, 'all_go_results')}, all_go_results={getattr(self, 'all_go_results', None)}")
            if hasattr(self, 'all_go_results') and self.all_go_results:
                print(f"Creating combined file with {len(self.all_go_results)} result sets")
                combined_go_results = pd.concat(self.all_go_results, ignore_index=True)
                combined_filename = "go_terms_with_genes_combined.csv"
                # combined_filename = get_final_path(combined_filename)
                combined_go_results.to_csv(combined_filename, index=False, encoding="utf-8")
                # Apply logic per cluster: for each cluster, select all significant if >5, otherwise top 5
                significant_go_results_list = []
                for cluster_val in combined_go_results['cluster'].unique():
                    cluster_data = combined_go_results[combined_go_results['cluster'] == cluster_val]
                    cluster_significant = cluster_data[cluster_data['sig'] != '-']
                    num_significant = cluster_significant.shape[0]
                    
                    if num_significant > 5:
                        # Select all significant ones for this cluster
                        significant_go_results_list.append(cluster_significant)
                    else:
                        # Select top 5 from all results for this cluster (sorted by p-value)
                        cluster_top5 = cluster_data.sort_values('p-value').head(5)
                        significant_go_results_list.append(cluster_top5)
                
                # Combine results from all clusters
                if significant_go_results_list:
                    significant_go_results = pd.concat(significant_go_results_list, ignore_index=True)
                else:
                    significant_go_results = pd.DataFrame()
                
                if significant_go_results.shape[0] > 0:
                    # Use base name from args.output and append _significant.csv
                    base_name = self.args.output.replace('.xlsx', '').replace('.xls', '').replace('.txt', '')
                    significant_filename = f"{base_name}_significant.csv"
                    # significant_filename = get_final_path(significant_filename)

                else:
                    print("No significant GO terms found (all have sig = '-')")
            else:
                print("No GO results to combine")
        if self.args.noprint == False:
            print("#"*80 + "\n")
        return outs

    def calcu_seurat_group(self,expfile,hgvc=False):
        """deal with seurat input matrix"""
        exps = read_csv(expfile)
        pre,suf,suf1 ="avg_logFC"," UMI counts/cell",""
        fid = "gene"
        pname = "p_val_adj"
        assert fid in exps.columns, 'No "gene" column. Wrong format? Seurat, Scanpy, Scran or Cellranger?'
        exps[fid] = exps[fid].str.replace("\.\d+","")
        cluster = "cluster"
        gcol = "gene" if hgvc == True else "ensemblID"
        ccol = "cellName"

        if self.args.target.lower() not in ["cancersea","cellmarker"]:
            print("Error target : -t, --target,(cellmarker,[cancersea])")
            sys.exit(0)

        if self.args.target.lower() == "cancersea":
            gcol = "gene" if hgvc == True else "ensemblID"
            ccol = "name"

        cnum = list(exps[cluster].unique())
        
        # Check if all cluster names are numerical
        all_numeric = all(str(x).isdigit() for x in cnum)
        
        if all_numeric:
            # Sort as numerical values
            cnum = sorted(cnum, key=lambda x: int(str(x)))
        else:
            # Sort by first character, then by the full string for secondary sorting
            cnum = sorted(cnum, key=lambda x: (str(x)[0], str(x)))
        abs_tag = True
        outs = []
        self.wb = self.wbgo = None
        if self.args.output:
            # Save to final folder
            original_output = self.args.output
            # self.args.output = get_final_path(self.args.output)
            if self.args.outfmt.lower() == "ms-excel":
                if not self.args.output.endswith(".xlsx") and (not self.args.output.endswith(".xls")):
                    self.args.output += ".xlsx"
                self.wb = ExcelWriter(self.args.output)
                self.wbgo = self.wb
            elif self.args.outfmt.lower() == "txt":
                self.wb = open(self.args.output,"w")
                if self.args.target == "cancersea":
                    self.wb.write("Cell Type\tZ-score\tNote\tCluster\n")
                else:
                    self.wb.write("Cell Type\tZ-score\tCluster\n")
                self.wbgo = open(self.args.output + ".go","w")
                self.wbgo.write('ids\tgene_num\tothergene_num\tp-value\tq-value\tsig\tname\tcluster\tgo_class\n')
            else:
                print("Error output format: -m, -outfmt,(ms-excel,[txt])")
                sys.exit(0)

        for i in cnum:
            cname = str(i)
            if self.args.cluster != "all":
                if self.args.cluster.find(",") > -1:
                    sets = self.args.cluster.split(",")
                    if cname not in sets:
                        continue
                else:
                    if cname != self.args.cluster:
                        continue
            o = " ".join(["#"*30,"Cluster",cname, "#"*30]) + "\n"
            if self.args.noprint == False:
                print(o)
            ltitle = pre
            ptitle = pname
            if ltitle not in exps.columns:
                print(ltitle,"column not in the input table!")
                sys.exit(0)
            newexps = exps[(exps[cluster] == i) & (exps[ltitle]>=self.args.foldchange) & (exps[ptitle] <= self.args.pvalue)]
            #newexps = exps[(exps[cluster] == i) & (abs(exps[ltitle])>=self.args.foldchange) & (exps[ptitle] <= self.args.pvalue)]
            #print(newexps)
            #print(newexps)

            h_values,colnames = self.get_cell_matrix(newexps,ltitle,fid,gcol,ccol,abs_tag)
            #print(colnames)
            #for x in newexps['gene'].unique():
            #    print(x)
            #exit()
            if self.args.output:
                h_values['Cluster'] = cname
                Annotator.to_output(h_values,self.wb,self.args.outfmt,cname,"Cell Type")

            #print(h_values)
            t,o_str,c,v,times = self.print_class(h_values,cname)
            outs.append([cname,t,c,v,times])
            if self.args.noprint == False:
                print(o_str)

            otherexps = exps[(exps[cluster] != i) & (exps[ltitle]>=self.args.foldchange) & (exps[ptitle] <= self.args.pvalue)]
            #otherexps = exps[(exps[cluster] != i) & (abs(exps[ltitle])>=self.args.foldchange) & (exps[ptitle] <= self.args.pvalue)]

            # if self.args.target.lower() == "cellmarker":
            # if 1==1:
            #     tfc,trownames,trownum,tcolnames,tcolnum = self.get_cell_gene_names(otherexps,self.cmarkers,fid,gcol,ccol,'other')
            #     if not trownames:continue
            #     other_gene_names = set(tcolnames)
            #     self.deal_with_badtype(cname,other_gene_names,colnames)
            # elif self.args.target.lower() == "cancersea":
            #     tfc,trownames,trownum,tcolnames,tcolnum = self.get_cell_gene_names(otherexps,self.smarkers,fid,gcol,ccol,'other')
            #     if not trownames:continue
            #     other_gene_names = set(tcolnames)
            #     self.deal_with_badtype(cname,other_gene_names,colnames)
            tfc,trownames,trownum,tcolnames,tcolnum = self.get_cell_gene_names(otherexps,self.cmarkers,fid,gcol,ccol,'other')
            # print(colnames)
            # input()
            self.deal_with_badtype(cname,set(tcolnames),colnames)
        if self.args.output:
            self.wb.close()
            self.wbgo.close()
            # Save combined GO results to single CSV file
            if hasattr(self, 'all_go_results') and self.all_go_results:
                combined_go_results = pd.concat(self.all_go_results, ignore_index=True)

                combined_filename = "go_terms_with_genes_combined.csv"
                # combined_filename = get_final_path(combined_filename)
                combined_go_results.to_csv(combined_filename, index=False, encoding="utf-8")
                # Apply logic per cluster: for each cluster, select all significant if >5, otherwise top 5
                significant_go_results_list = []
                for cluster_val in combined_go_results['cluster'].unique():
                    cluster_data = combined_go_results[combined_go_results['cluster'] == cluster_val]
                    cluster_significant = cluster_data[cluster_data['sig'] != '-']
                    num_significant = cluster_significant.shape[0]
                    
                    if num_significant > 5:
                        # Select all significant ones for this cluster
                        significant_go_results_list.append(cluster_significant)
                    else:
                        # Select top 5 from all results for this cluster (sorted by p-value)
                        cluster_top5 = cluster_data.sort_values('p-value').head(5)
                        significant_go_results_list.append(cluster_top5)
                
                # Combine results from all clusters
                if significant_go_results_list:
                    significant_go_results = pd.concat(significant_go_results_list, ignore_index=True)
                else:
                    significant_go_results = pd.DataFrame()
                
                if significant_go_results.shape[0] > 0:
                    # Use base name from args.output and append _significant.csv
                    base_name = self.args.output.replace('.xlsx', '').replace('.xls', '').replace('.txt', '')
                    significant_filename = f"{base_name}_significant.csv"
                    # significant_filename = get_final_path(significant_filename)
                    significant_go_results.to_csv(significant_filename, index=False, encoding="utf-8")
                else:
                    print("No significant GO terms found (all have sig = '-')")
        if self.args.noprint == False:
            print("#"*80 + "\n")
        return outs

    def calcu_scanpy_group(self,expfile,hgvc=False):
        """deal with scanpy input matrix"""
        exps = read_csv(expfile,index_col=0)
        cnum = set()
        pname = "p"
        pre = "l"
        fid = "n"
        for c in exps.columns:
            k,v = c.split("_")
            cnum.add(k)
            if v.startswith("p"):
                pname = v
            elif v.startswith("n"):
                rfid = v
            elif v.startswith("l"):
                pre = v
        
        #pre,suf,suf1 ="avg_logFC"," UMI counts/cell",""
        #fid = "gene"
        #pname = "p_val_adj"
        #assert fid in exps.columns, 'No "gene" column. Wrong format? Scanpy, Seurat or Cellranger?'
        #exps[fid] = exps[fid].str.replace("\.\d+","")
        #cluster = "cluster"

        ###MarkerBase
        gcol = "gene" if hgvc == True else "ensemblID"
        ccol = "cellName"

        if self.args.target.lower() not in ["cancersea","cellmarker"]:
            print("Error target : -t, --target,(cellmarker,[cancersea])")
            sys.exit(0)

        if self.args.target.lower() == "cancersea":
            gcol = "gene" if hgvc == True else "ensemblID"
            ccol = "name"

        # Apply same cluster name ordering as seurat
        cnum = list(cnum)
        
        # Check if all cluster names are numerical
        all_numeric = all(str(x).isdigit() for x in cnum)
        
        if all_numeric:
            # Sort as numerical values
            cnum = sorted(cnum, key=lambda x: int(str(x)))
        else:
            # Sort by first character, then by the full string for secondary sorting
            cnum = sorted(cnum, key=lambda x: (str(x)[0], str(x)))
        
        abs_tag = True
        outs = []
        self.wb = self.wbgo = None
        if self.args.output:
            # Save to final folder
            original_output = self.args.output
            # self.args.output = get_final_path(self.args.output)
            if self.args.outfmt.lower() == "ms-excel":
                if not self.args.output.endswith(".xlsx") and (not self.args.output.endswith(".xls")):
                    self.args.output += ".xlsx"
                self.wb = ExcelWriter(self.args.output)
                self.wbgo = self.wb
            elif self.args.outfmt.lower() == "txt":
                self.wb = open(self.args.output,"w")
                if self.args.target == "cancersea":
                    self.wb.write("Cell Type\tZ-score\tNote\tCluster\n")
                else:
                    self.wb.write("Cell Type\tZ-score\tCluster\n")
                self.wbgo = open(self.args.output + ".go","w")
                self.wbgo.write('ids\tgene_num\tothergene_num\tp-value\tq-value\tsig\tname\tcluster\tgo_class\n')
            else:
                print("Error output format: -m, -outfmt,(ms-excel,[txt])")
                sys.exit(0)

        for i in cnum:
            cname = str(i)
            if self.args.cluster != "all":
                if self.args.cluster.find(",") > -1:
                    sets = self.args.cluster.split(",")
                    if cname not in sets:
                        continue
                else:
                    if cname != self.args.cluster:
                        continue
            o = " ".join(["#"*30,"Cluster",cname, "#"*30]) + "\n"
            if self.args.noprint == False:
                print(o)
            ltitle = cname + "_" + pre
            fid = cname + "_" + rfid
            ptitle = cname + "_" + pname
            if ltitle not in exps.columns:
                print(ltitle,"column not in the input table!")
                sys.exit(0)
            newexps = exps[[fid,ltitle,ptitle]][(exps[ltitle]>=self.args.foldchange) & (exps[ptitle] <= self.args.pvalue)]
            #newexps = exps[(exps[cluster] == i) & (abs(exps[ltitle])>=self.args.foldchange) & (exps[ptitle] <= self.args.pvalue)]
            #print(newexps)
            #print(newexps)

            h_values,colnames = self.get_cell_matrix(newexps,ltitle,fid,gcol,ccol,abs_tag)
            print("Cluster " + cname + " Gene number:",newexps[fid].unique().shape[0])
            #print(colnames)
            #for x in newexps[fid].unique():
            #    print(x)
            #exit()
            if self.args.output:
                h_values['Cluster'] = cname
                Annotator.to_output(h_values,self.wb,self.args.outfmt,cname,"Cell Type")

            #print(h_values)
            #exit()
            t,o_str,c,v,times = self.print_class(h_values,cname)
            outs.append([cname,t,c,v,times])
            if self.args.noprint == False:
                print(o_str)

            otherexps = None
            ofid = 'o_n'
            oltitle = 'o_l'
            optitle = 'o_p'
            for j in cnum:
                oname = str(j)
                if oname == cname:continue
                tltitle = oname + "_" + pre
                tfid = oname + "_" + rfid
                tptitle = oname + "_" + pname
                tempexps = exps[[tfid,tltitle,tptitle]][(exps[tltitle]>=self.args.foldchange) & (exps[tptitle] <= self.args.pvalue)]
                tempexps.columns = [ofid,oltitle,optitle]
                if otherexps is None:
                    otherexps = tempexps
                else:
                    otherexps = pd.concat([otherexps,tempexps])
            #otherexps = exps[(exps[cluster] != i) & (abs(exps[ltitle])>=self.args.foldchange) & (exps[ptitle] <= self.args.pvalue)]
            #print(otherexps)
            #exit()

            if self.args.target.lower() == "cellmarker":
                tfc,trownames,trownum,tcolnames,tcolnum = self.get_cell_gene_names(otherexps,self.cmarkers,ofid,gcol,ccol,'other')
                if not trownames:continue
                other_gene_names = set(tcolnames)
                self.deal_with_badtype(cname,other_gene_names,colnames)
            elif self.args.target.lower() == "cancersea":
                tfc,trownames,trownum,tcolnames,tcolnum = self.get_cell_gene_names(otherexps,self.smarkers,ofid,gcol,ccol,'other')
                if not trownames:continue
                other_gene_names = set(tcolnames)
                self.deal_with_badtype(cname,other_gene_names,colnames)
            print("Other Gene number:",len(other_gene_names))
        if self.args.output:
            self.wb.close()
            self.wbgo.close()
            # Save combined GO results to single CSV file
            if hasattr(self, 'all_go_results') and self.all_go_results:
                combined_go_results = pd.concat(self.all_go_results, ignore_index=True)
                # Convert cluster column to numeric for proper numerical sorting
                combined_go_results['cluster'] = pd.to_numeric(combined_go_results['cluster'], errors='coerce')
                combined_go_results = combined_go_results.sort_values(['cluster', 'p-value'])
                combined_filename = "go_terms_with_genes_combined.csv"
                # combined_filename = get_final_path(combined_filename)
                combined_go_results.to_csv(combined_filename, index=False, encoding="utf-8")
                # Apply logic per cluster: for each cluster, select all significant if >5, otherwise top 5
                significant_go_results_list = []
                for cluster_val in combined_go_results['cluster'].unique():
                    cluster_data = combined_go_results[combined_go_results['cluster'] == cluster_val]
                    cluster_significant = cluster_data[cluster_data['sig'] != '-']
                    num_significant = cluster_significant.shape[0]
                    
                    if num_significant > 5:
                        # Select all significant ones for this cluster
                        significant_go_results_list.append(cluster_significant)
                    else:
                        # Select top 5 from all results for this cluster (sorted by p-value)
                        cluster_top5 = cluster_data.sort_values('p-value').head(5)
                        significant_go_results_list.append(cluster_top5)
                
                # Combine results from all clusters
                if significant_go_results_list:
                    significant_go_results = pd.concat(significant_go_results_list, ignore_index=True)
                else:
                    significant_go_results = pd.DataFrame()
                
                if significant_go_results.shape[0] > 0:
                    # Use base name from args.output and append _significant.csv
                    base_name = self.args.output.replace('.xlsx', '').replace('.xls', '').replace('.txt', '')
                    significant_filename = f"{base_name}_significant.csv"
                    # significant_filename = get_final_path(significant_filename)
                    significant_go_results.to_csv(significant_filename, index=False, encoding="utf-8")
                else:
                    print("No significant GO terms found (all have sig = '-')")
        if self.args.noprint == False:
            print("#"*80 + "\n")
        return outs

    def calcu_scran_group(self,expfile,hgvc=False):
        """deal with scran input matrix"""
        exps = read_csv(expfile)
        if exps.columns[0] == "Unnamed: 0":
            exps.rename(columns={"Unnamed: 0":"gene"},inplace=True)
        cnum = set()
        pname = "p.value"
        pre = "LFC"
        fid = "gene"
        for c in exps.columns:
            if c.startswith("LFC"):
                k,v = c.split("_")
                cnum.add(v)
        
        #pre,suf,suf1 ="avg_logFC"," UMI counts/cell",""
        #fid = "gene"
        #pname = "p_val_adj"
        #assert fid in exps.columns, 'No "gene" column. Wrong format? Scanpy, Seurat or Cellranger?'
        #exps[fid] = exps[fid].str.replace("\.\d+","")
        #cluster = "cluster"

        ###MarkerBase
        gcol = "gene" if hgvc == True else "ensemblID"
        ccol = "cellName"

        if self.args.target.lower() not in ["cancersea","cellmarker"]:
            print("Error target : -t, --target,(cellmarker,[cancersea])")
            sys.exit(0)

        if self.args.target.lower() == "cancersea":
            gcol = "gene" if hgvc == True else "ensemblID"
            ccol = "name"

        #cnum = list(exps[cluster].unique())
        abs_tag = True
        outs = []
        self.wb = self.wbgo = None
        if self.args.output:
            # Save to final folder
            original_output = self.args.output
            # self.args.output = get_final_path(self.args.output)
            if self.args.outfmt.lower() == "ms-excel":
                if not self.args.output.endswith(".xlsx") and (not self.args.output.endswith(".xls")):
                    self.args.output += ".xlsx"
                self.wb = ExcelWriter(self.args.output)
                self.wbgo = self.wb
            elif self.args.outfmt.lower() == "txt":
                self.wb = open(self.args.output,"w")
                if self.args.target == "cancersea":
                    self.wb.write("Cell Type\tZ-score\tNote\tCluster\n")
                else:
                    self.wb.write("Cell Type\tZ-score\tCluster\n")
                self.wbgo = open(self.args.output + ".go","w")
                self.wbgo.write('ids\tgene_num\tothergene_num\tp-value\tq-value\tsig\tname\tcluster\tgo_class\n')
            else:
                print("Error output format: -m, -outfmt,(ms-excel,[txt])")
                sys.exit(0)

        for i in list(sorted(cnum)):
            cname = str(i)
            if self.args.cluster != "all":
                if self.args.cluster.find(",") > -1:
                    sets = self.args.cluster.split(",")
                    if cname not in sets:
                        continue
                else:
                    if cname != self.args.cluster:
                        continue
            o = " ".join(["#"*30,"Cluster",cname, "#"*30]) + "\n"
            if self.args.noprint == False:
                print(o)
            ltitle = pre + "_" + cname
            ptitle = pname + "_" + cname
            if ltitle not in exps.columns:
                print(ltitle,"column not in the input table!")
                sys.exit(0)
            newexps = exps[[fid,ltitle,ptitle]][(exps[ltitle]>=self.args.foldchange) & (exps[ptitle] <= self.args.pvalue)]
            #newexps = exps[(exps[cluster] == i) & (abs(exps[ltitle])>=self.args.foldchange) & (exps[ptitle] <= self.args.pvalue)]
            #print(newexps)
            #print(newexps)
            #continue

            h_values,colnames = self.get_cell_matrix(newexps,ltitle,fid,gcol,ccol,abs_tag)
            #print(colnames)
            #for x in newexps[fid].unique():
            #    print(x)
            #exit()
            if self.args.output:
                h_values['Cluster'] = cname
                Annotator.to_output(h_values,self.wb,self.args.outfmt,cname,"Cell Type")

            #print(h_values)
            #exit()
            t,o_str,c,v,times = self.print_class(h_values,cname)
            outs.append([cname,t,c,v,times])
            if self.args.noprint == False:
                print(o_str)

            otherexps = None
            ofid = 'o_n'
            oltitle = 'o_l'
            optitle = 'o_p'
            for j in list(sorted(cnum)):
                oname = str(j)
                if oname == cname:continue
                tltitle = pre + "_" + oname
                tfid = fid 
                tptitle = pname + "_" + oname
                tempexps = exps[[tfid,tltitle,tptitle]][(exps[tltitle]>=self.args.foldchange) & (exps[tptitle] <= self.args.pvalue)]
                tempexps.columns = [ofid,oltitle,optitle]
                if otherexps is None:
                    otherexps = tempexps
                else:
                    otherexps = pd.concat([otherexps,tempexps])
            #otherexps = exps[(exps[cluster] != i) & (abs(exps[ltitle])>=self.args.foldchange) & (exps[ptitle] <= self.args.pvalue)]
            #print(otherexps)
            #exit()

            if self.args.target.lower() == "cellmarker":
                tfc,trownames,trownum,tcolnames,tcolnum = self.get_cell_gene_names(otherexps,self.cmarkers,ofid,gcol,ccol,'other')
                if not trownames:continue
                other_gene_names = set(tcolnames)
                self.deal_with_badtype(cname,other_gene_names,colnames)
            elif self.args.target.lower() == "cancersea":
                tfc,trownames,trownum,tcolnames,tcolnum = self.get_cell_gene_names(otherexps,self.smarkers,ofid,gcol,ccol,'other')
                if not trownames:continue
                other_gene_names = set(tcolnames)
                self.deal_with_badtype(cname,other_gene_names,colnames)
            print("Other Gene number:",len(other_gene_names))
        if self.args.output:
            self.wb.close()
            self.wbgo.close()
            # Save combined GO results to single CSV file
            if hasattr(self, 'all_go_results') and self.all_go_results:
                combined_go_results = pd.concat(self.all_go_results, ignore_index=True)
                # Convert cluster column to numeric for proper numerical sorting
                combined_go_results['cluster'] = pd.to_numeric(combined_go_results['cluster'], errors='coerce')
                combined_go_results = combined_go_results.sort_values(['cluster', 'p-value'])
                combined_filename = "go_terms_with_genes_combined.csv"
                # combined_filename = get_final_path(combined_filename)
                combined_go_results.to_csv(combined_filename, index=False, encoding="utf-8")
                # Apply logic per cluster: for each cluster, select all significant if >5, otherwise top 5
                significant_go_results_list = []
                for cluster_val in combined_go_results['cluster'].unique():
                    cluster_data = combined_go_results[combined_go_results['cluster'] == cluster_val]
                    cluster_significant = cluster_data[cluster_data['sig'] != '-']
                    num_significant = cluster_significant.shape[0]
                    
                    if num_significant > 5:
                        # Select all significant ones for this cluster
                        significant_go_results_list.append(cluster_significant)
                    else:
                        # Select top 5 from all results for this cluster (sorted by p-value)
                        cluster_top5 = cluster_data.sort_values('p-value').head(5)
                        significant_go_results_list.append(cluster_top5)
                
                # Combine results from all clusters
                if significant_go_results_list:
                    significant_go_results = pd.concat(significant_go_results_list, ignore_index=True)
                else:
                    significant_go_results = pd.DataFrame()
                
                if significant_go_results.shape[0] > 0:
                    # Use base name from args.output and append _significant.csv
                    base_name = self.args.output.replace('.xlsx', '').replace('.xls', '').replace('.txt', '')
                    significant_filename = f"{base_name}_significant.csv"
                    # significant_filename = get_final_path(significant_filename)
                    significant_go_results.to_csv(significant_filename, index=False, encoding="utf-8")
                else:
                    print("No significant GO terms found (all have sig = '-')")
        if self.args.noprint == False:
            print("#"*80 + "\n")
        return outs


    def get_exp_matrix_loop(self,exps,ltitle,fid,colnames,rownames,cell_matrix,usertag,abs_tag = True):
        """format the cell_deg_matrix and calculate the zscore of certain cell types."""

        #filter gene expressed matrix according to the markers
        gene_exps = exps.loc[:,[fid,ltitle]][exps[fid].isin(colnames)]

        gene_matrix = mat(gene_exps.sort_values(fid)[ltitle]).T
        gene_matrix = gene_matrix * np.mean(gene_matrix) ### / np.min(gene_matrix))

        if gene_matrix.shape[0] != cell_matrix.shape[1]:
            #print(gene_matrix.shape,cell_matrix.shape)
            #print(len(gene_exps[fid].unique()))
            #print(gene_matrix)
            #print(cell_matrix)
            print("Error for inconsistent gene numbers, please check your expression csv for '" + fid + "'")
            return None
        
        nonzero = np.matrix(np.count_nonzero(cell_matrix,axis=1)).T
        #gene_matrix = np.ones_like(gene_matrix)
        cell_deg_matrix = cell_matrix * gene_matrix

        #print("cell",cell_matrix)
        #print("gene",gene_matrix)
        #print(colnames)
        #print(rownames)

        #print(rownames)
        #exit()
        #print(gene_matrix)
        #print(cell_deg_matrix)
        #print(type(rownames))
        #a1 = "Natural killer T (NKT) cell"
        #b1 = "T cell"
        #a1 = "Macrophage"
        #b1 = "Monocyte"

        #a1 = "Mesenchymal stem cell"
        #b1 = "Fibroblast"
        #mar = cell_matrix[np.array(rownames) == a1]
        #mon = cell_matrix[np.array(rownames) == b1]
        #marz = nonzero[np.array(rownames) == a1]
        #monz = nonzero[np.array(rownames) == b1]
        #print(len(mar[np.nonzero(mar)]),len(mon[np.nonzero(mon)]))

        #print(log2(marz),log2(monz))

        #print(mar)
        #print(mon)
        #print(marz,monz)
        #print(cell_matrix,cell_matrix.shape,gene_matrix.shape)

        #print(np.std(cell_matrix,axis=1))
        #print(cell_matrix.shape,cell_deg_matrix.shape)
        wstd = np.matrix(np.std(cell_matrix,axis=1)).T
        #print(wstd.shape,wstd,nonzero)
        if usertag:
            cell_deg_matrix = np.matrix(np.array(cell_deg_matrix))
        else:
            if (wstd.shape == np.ones_like(wstd)).all:
                wstd = [[1]]
            if (nonzero == np.ones_like(nonzero)).all:
                cell_deg_matrix = np.matrix(np.array(cell_deg_matrix) * np.array(wstd))
            else:
                cell_deg_matrix = np.matrix(np.array(cell_deg_matrix) * np.array(log2(nonzero)) * np.array(wstd))

        out = DataFrame({"Z-score":cell_deg_matrix.A1},index=rownames)
        out.sort_values(['Z-score'],inplace=True,ascending=False)
        #out.to_csv("wei.sco",sep="\t")
        #print(cell_deg_matrix,wstd,log2(nonzero))

        if abs_tag:
            out['Z-score'] = abs(out['Z-score'])
        else:
            out = out[out['Z-score'] > 0]

        #print(out)
        if (out.shape[0] > 1):
            out['Z-score'] = (out['Z-score'] - mean(out['Z-score']))/std(out['Z-score'],ddof=1)
        #print(out)

        return out


    def get_cell_gene_names(self,exps,markers,fid,gcol,ccol,tag):
        """find expressed markers according to the markers and expressed matrix."""
        #print(fid)
        whole_gsets = set(exps[fid])
        if self.args.target.lower() == "cancersea":
            #whole_fil = markers['EnsembleID'].isin(whole_gsets)
            markers['weight'] = 1
            if self.args.Gensymbol == True:
                markers[gcol] = markers['GeneName']
            else:
                markers[gcol] = markers['EnsembleID']
        whole_fil = markers[gcol].isin(whole_gsets)

        fc = markers[[ccol,gcol,'weight']][whole_fil]
        #print(whole_gsets,exps[fid],exps)
        #print(markers[gcol])
        #print(fc)
        #print(list(fc['cellName'].unique()),ccol,gcol)
        #print(whole_fil.unique())
        #exit()
        #print(markers,markers.columns)
        #print(fc)
        if fc.shape[0] == 0:
            if tag != "other":
                print("!WARNING3:Zero marker sets found, type:" + tag)
                print("!WARNING3:Change the threshold or tissue name and try again?")
                print("!WARNING3:EnsemblID or GeneID,try '-E' command?")
            return fc,None,None,whole_gsets,None
        #print("helll")
        fc.columns = [ccol,gcol,'c']
        fc.set_index([ccol,gcol])
        newfc = fc.groupby([ccol,gcol]).sum()
        #print(newfc)
        #if newfc.shape[0] <1:
        #    print(newfc.shape)
        #    print(fc)
        #    print(exps)
        names = newfc.index
        #print(names)
        #print(names)
        newfc['c1'] = names
        newfc[gcol] = newfc['c1'].apply(lambda x:x[1])
        newfc[ccol] = newfc['c1'].apply(lambda x:x[0])
        newfc.drop(['c1'],inplace=True,axis=1)
        newfc.reset_index(drop=True,inplace=True)
        #print(newfc)
        #exit()
        newfc['c'] = log2(newfc['c'] + 0.05) # * np.min(newfc['c'])
        fc = newfc
        #print("hello")
        #newfc.to_csv("wei.cls",sep="\t")
        #exit()
        #print(fc['c'][fc['c'] != 0])

        rownames = sorted(set(fc[ccol].unique()))
        rownum = len(rownames)
        colnames = sorted(set(fc[gcol].unique()))
        colnum = len(colnames)
        #print(fc.shape,fc)
        return fc,rownames,rownum,colnames,colnum

    def get_user_cell_gene_names(self,exps,fid,gcol,ccol,tag):
        """find expressed markers according to the user markers and expressed matrix."""
        #print(self.usermarkers)
        self.usermarkers.columns = [ccol,gcol,'weight']
        whole_gsets = set(exps[fid])
        whole_fil = self.usermarkers[gcol].isin(whole_gsets)

        fc = self.usermarkers[[ccol,gcol,'weight']][whole_fil]
        if fc.shape[0] == 0:
            if tag != "other":
                print("!WARNING3:Zero marker sets found, type:" + tag)
                print("!WARNING3:Change the threshold or tissue name and try again?")
                print("!WARNING3:EnsemblID or GeneID,try '-E' command?")
            return fc,None,None,whole_gsets,None
        fc.columns = [ccol,gcol,'c']
        fc.set_index([ccol,gcol])
        #print("FC",fc)
        #print("ENSG00000105369" in whole_gsets)

        newfc = fc.groupby([ccol,gcol]).sum()
        #if newfc.shape[0] <1:
        #    print(newfc.shape)
        #    print(fc)
        #    print(exps)
        names = newfc.index
        #print(names)
        #print(names)
        newfc['c1'] = names
        newfc[gcol] = newfc['c1'].apply(lambda x:x[1])
        newfc[ccol] = newfc['c1'].apply(lambda x:x[0])
        newfc.drop(['c1'],inplace=True,axis=1)
        newfc.reset_index(drop=True,inplace=True)
        #print(newfc)
        #exit()
        newfc['c'] = log2(newfc['c'] + 0.05) # * np.min(newfc['c'])
        fc = newfc
        #print("hello")
        #newfc.to_csv("wei.cls",sep="\t")
        #exit()


        rownames = sorted(set(self.usermarkers[ccol].unique()))
        rownum = len(rownames)
        colnames = sorted(set(fc[gcol].unique()))
        colnum = len(colnames)

        #print(fc,rownames,colnames)
        #print(fc.shape)
        return fc,rownames,rownum,colnames,colnum

    def get_cell_matrix(self,exps,ltitle,fid,gcol,ccol,abs_tag):
        """combine cell matrix with weight-matrix"""
        cell_value = None
        colnames = None
        if not self.args.norefdb:
            cell_value,colnames =self.get_cell_matrix_detail(exps,ltitle,fid,gcol,ccol,False,abs_tag)
        if self.args.MarkerDB != None:
            if self.args.norefdb:
                cell_value,colnames =self.get_cell_matrix_detail(exps,ltitle,fid,gcol,ccol,True,abs_tag)
            else:
                cell_value,colnames =self.get_cell_matrix_detail(exps,ltitle,fid,gcol,ccol,False,abs_tag)
                user_value,user_colnames =self.get_cell_matrix_detail(exps,ltitle,fid,gcol,ccol,True,abs_tag)
                #print("C",cell_value)
                #print("U",user_value)
                if cell_value is None:
                    if user_value is None:
                        if colnames is None:
                            return DataFrame(),None
                        else:
                            return DataFrame(),set(colnames)
                    else:
                        cell_value = user_value
                        colnames = user_colnames
                        cell_value = cell_value.join(user_value,how="outer",lsuffix="cm",rsuffix="ur")
                        cell_value[cell_value.isna()] = 0
                        colnames = colnames | user_colnames
                else:
                    if user_value is None:
                        user_value = cell_value
                        user_colnames = colnames
                    cell_value = cell_value.join(user_value,how="outer",lsuffix="cm",rsuffix="ur")
                    cell_value[cell_value.isna()] = 0
                    colnames = colnames | user_colnames

                #else:
                #    cell_value = user_value
                #    colnames = user_colnames
        #database weight-matrix
        wm = [1]
        if self.args.MarkerDB != None:
            if self.args.norefdb:
                wm = [1]
            else:
                wm =[0.1,0.9]
        weight_matrix = mat(wm).T

        if colnames is None:
            return DataFrame(),None

        if cell_value is None:
            return DataFrame(),set(colnames)

        #print(cell_value)
        last_value = array(cell_value) * weight_matrix
        result = DataFrame({"Cell Type":cell_value.index,"Z-score":last_value.A1})
        result = result.sort_values(by="Z-score",ascending = False)
        #if self.args.target == "cancersea":
        #    result['note'] = result['Cell Type'].apply(lambda x: self.snames[x])
        return result,set(colnames)

    def get_cell_matrix_detail(self,exps,ltitle,fid,gcol,ccol,usertag,abs_tag):
        """calculate the cell type scores"""
        fc,rownames,rownum,colnames,colnum = None,None,None,None,None
        #print(self.cmarkers)
        if self.args.target == "cellmarker":
            markers = self.cmarkers
        elif self.args.target == "cancersea":
            markers = self.smarkers

        #print(markers.columns)

        if usertag:
            fc,rownames,rownum,colnames,colnum = self.get_user_cell_gene_names(exps,fid,gcol,ccol,"user_marker")
            #print("FC",fc)
            #fc['c'] = 1
        else:
            fc,rownames,rownum,colnames,colnum = self.get_cell_gene_names(exps,markers,fid,gcol,ccol,'marker')
            #print(colnames)
        #print(colnames)

        if not colnames:
            return None,None
        if fc.shape[0] == 0:
            return None,set(colnames)

        exps = exps[exps[fid].isin(colnames)]

        rowdic = dict(zip(rownames,range(rownum)))
        coldic = dict(zip(colnames,range(colnum)))
        fc_cell = fc[ccol].map(lambda x:rowdic[x])
        fc_gene = fc[gcol].map(lambda x:coldic[x])

        newdf = DataFrame({ccol:fc_cell,gcol:fc_gene,"c":fc['c']})
        cell_coo_matrix = coo_matrix((newdf['c'],(newdf[ccol],newdf[gcol])),shape=(rownum,colnum))
        cell_matrix = cell_coo_matrix.toarray()

        #print(newdf)
        #print(rownames)
        #print(colnames)
        #print(cell_matrix)

        if self.args.noprint == False:
            if usertag:
                print("User Cell Num:",rownum)
                print("User Gene Num:",colnum)
                print("User Not Zero:",cell_coo_matrix.count_nonzero())
            else:
                print("Cell Num:",rownum)
                print("Gene Num:",colnum)
                print("Not Zero:",cell_coo_matrix.count_nonzero())
        cell_values = self.get_exp_matrix_loop(exps,ltitle,fid,colnames,rownames,cell_matrix,usertag,abs_tag)
        return cell_values,set(colnames)


    def read_user_markers(self,colname):
        """usermarker db preparation"""
        if self.args.MarkerDB != None:
            if not os.path.exists(self.args.MarkerDB):
                print("User marker database does not exists!",self.args.MarkerDB)
                sys.exit(0)
            self.usermarkers = read_csv(self.args.MarkerDB,sep="\t",header=None)
            self.usermarkers.columns=['cellName',colname]
            #self.hgncs_ensem = dict(zip(self.ensem_hgncs.values(),self.ensem_hgncs.keys()))
            if colname == "ensemblID":
                self.usermarkers[colname] = self.usermarkers[colname].map(lambda x:self.hgncs_ensem[x] if x in self.hgncs_ensem else x)
            self.usermarkers['weight'] = 1
            if self.args.noprint == False:
                print("User cells:", len(self.usermarkers['cellName'].unique()))
                print("User genes:", len(self.usermarkers[colname].unique()))

    def load_pickle_module(self,db):
        """read whole database"""
        handler = gzip.open(db,"rb")
        self.gos = load(handler)
        self.human_gofs = load(handler)
        self.mouse_gofs = load(handler)
        self.cmarkers = load(handler)
        self.smarkers = load(handler)
        self.snames = load(handler)
        self.ensem_hgncs = load(handler)
        self.ensem_mouse = load(handler)
        self.hgncs_ensem = dict(zip(self.ensem_hgncs.values(),self.ensem_hgncs.keys()))
        fil = []
        #fil = ['Cancer stem cell', 'Cancer cell']
        #print(self.cmarkers)
        #exit()
        self.cmarkers = self.cmarkers[~self.cmarkers['cellName'].isin(fil)]

        #if self.args.noprint == False:
        print("Version V1.1 [2020/07/03]")
        print("DB load:",len(self.gos),len(self.human_gofs),len(self.mouse_gofs),len(self.cmarkers),len(self.ensem_hgncs))

    def read_tissues_species(self,tissue="All",species="Human",celltype="normal"):
        """read markers according to certain tissue and certain species"""
        species = species.lower().capitalize()
        ct = celltype.lower().capitalize()
        if tissue != "All":
            if tissue.find(",") > -1:
                temptis = set(tissue.split(","))
                self.cmarkers = self.cmarkers[self.cmarkers['tissueType'].isin(temptis)]
            else:
                self.cmarkers = self.cmarkers[self.cmarkers['tissueType'].isin([tissue])]
        if ct == "Normal":
            self.cmarkers = self.cmarkers[self.cmarkers['cellType']=="Normal cell"]
        elif ct == "Cancer":
            self.cmarkers = self.cmarkers[self.cmarkers['cellType']=="Cancer cell"]
        else:
            print("Illegal celltype. Please use \"[Normal] or [Cancer] instead.")
            exit(0)


        #self.cmarkers = self.cmarkers[self.cmarkers['cellName']!="Mesenchymal stem cell"]
        print("load markers:",len(self.cmarkers))
        self.cmarkers = self.cmarkers[self.cmarkers['speciesType'].isin([species])]
        #print(self.cmarkers)

    def get_list_tissue(self,species):
        """print tissue names"""
        species = species.lower().capitalize()
        cmarkers = self.cmarkers[self.cmarkers['speciesType'].isin([species])]
        names = list(sorted(cmarkers['tissueType'].unique()))
        print("#" * 120)
        print("-" * 120)
        print("{0:s}{1:<10s}{2:>5s}{3:<10d}".format("Species:",species,"Num:",len(names)))
        print("-" * 120)
        for i in range(0,len(names)-2,3):
            if len(names) < i + 1:
                s = "{0:3d}: {1:<40s}".format(i+1,names[i])
            elif len(names) < i + 2:
                s = "{0:3d}: {1:<35s}{2:3d}: {3:<35s}".format(i+1,names[i],i+2,names[i+1])
            else:
                s = "{0:3d}: {1:<35s}{2:3d}: {3:<35s}{4:3d}: {5:<35s}".format(i+1,names[i],i+2,names[i+1],i+3,names[i+2])
            print(s)
        print("#" * 120)




    def run_detail_cmd(self):
        """main command"""
        #self.check_db()
        if not os.path.exists(self.args.input):
            tempname = "./" + self.args.input
            if not os.path.exists(tempname):
                print(tempname)
                print("Input file does not exists!",self.args.input)
                sys.exit(0)
        print(self.args)
        if self.args.source.lower() == "cellranger":
            self.load_pickle_module(self.args.db)
            if self.args.species == "Mouse":
                self.ensem_hgncs = self.ensem_mouse
                self.human_gofs = self.mouse_gofs
            self.read_tissues_species(self.args.tissue,self.args.species,self.args.celltype)
            if self.args.Gensymbol:
                self.read_user_markers('Gene ID')
            else:
                self.read_user_markers('ensemblID')
            outs = self.calcu_cellranger_group(self.args.input,self.args.Gensymbol)
            return outs
        elif args.source.lower() == "seurat":
            self.load_pickle_module(self.args.db)
            if self.args.species == "Mouse":
                self.ensem_hgncs = self.ensem_mouse
                self.human_gofs = self.mouse_gofs
            self.read_tissues_species(self.args.tissue,self.args.species,self.args.celltype)
            if self.args.Gensymbol:
                self.read_user_markers('gene')
            else:
                self.read_user_markers('ensemblID')
            outs = self.calcu_seurat_group(self.args.input,self.args.Gensymbol)
            return outs
        elif args.source.lower() == "scanpy":
            self.load_pickle_module(self.args.db)
            if self.args.species == "Mouse":
                self.ensem_hgncs = self.ensem_mouse
                self.human_gofs = self.mouse_gofs
            self.read_tissues_species(self.args.tissue,self.args.species,self.args.celltype)
            if self.args.Gensymbol:
                self.read_user_markers('gene')
            else:
                self.read_user_markers('ensemblID')
            outs = self.calcu_scanpy_group(self.args.input,self.args.Gensymbol)
            return outs
        elif args.source.lower() == "scran":
            self.load_pickle_module(self.args.db)
            if self.args.species == "Mouse":
                self.ensem_hgncs = self.ensem_mouse
                self.human_gofs = self.mouse_gofs
            self.read_tissues_species(self.args.tissue,self.args.species,self.args.celltype)
            if self.args.Gensymbol:
                self.read_user_markers('gene')
            else:
                self.read_user_markers('ensemblID')
            outs = self.calcu_scran_group(self.args.input,self.args.Gensymbol)
            return outs
            pass



class Process(object):
    def __init__(self):
        pass
    def get_parser(self):
        desc = """Program: SCSA
  Version: 1.0
  Email  : <yhcao@ibms.pumc.edu.cn>
        """
        #usage = "%(prog)s"
        parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=desc)
                                         #usage=usage)

        parser.add_argument('-i', '--input', required = True, help="Input file for marker annotation(Only CSV format supported).")
        parser.add_argument('-o', '--output', help="Output file for marker annotation.")
        parser.add_argument('-d', '--db', default = "whole.db",help="Database for annotation. (whole.db)")
        parser.add_argument('-s', '--source', default = "cellranger",help="Source of marker genes. (cellranger,[seurat],[scanpy],[scran])")
        parser.add_argument('-c', '--cluster', default = "all",help="Only deal with one cluster of marker genes. (all,[1],[1,2,3],[...])")
        parser.add_argument('-M', '--MarkerDB', help='User-defined marker database in table format with two columns.First column as Cellname, Second refers to Genename.')
        parser.add_argument('-f',"--foldchange",default = 2,help="Fold change threshold for marker filtering. (2.0)")
        parser.add_argument('-p',"--pvalue",default = 0.05,help="P-value threshold for marker filtering. (0.05)")
        parser.add_argument('-w',"--weight",default = 100,help="Weight threshold for marker filtering from cellranger v1.0 results. (100)")
        parser.add_argument('-g',"--species",default = 'Human',help="Species for annotation. Only used for cellmarker database. ('Human',['Mouse'])")
        parser.add_argument('-k',"--tissue",default = 'All',help="Tissue for annotation. Only used for cellmarker database. Multiple tissues should be seperated by commas. Run '-l' option to see all tissues. ('All',['Bone marrow'],['Bone marrow,Brain,Blood'][...])")
        parser.add_argument('-m', '--outfmt', default = "ms-excel", help="Output file format for marker annotation. (ms-excel,[txt])")
        parser.add_argument('-T',"--celltype",default = "normal",help="Cell type for annotation. (normal,[cancer])")
        parser.add_argument('-t', '--target', default = "cellmarker",help="Target to annotation class in Database. (cellmarker,[cancersea])")
        parser.add_argument('-E',"--Gensymbol",action = "store_true",default=False,help="Using gene symbol ID instead of ensembl ID in input file for calculation.")
        parser.add_argument('-N',"--norefdb",action = "store_true",default=False,help="Only using user-defined marker database for annotation.")
        parser.add_argument('-b',"--noprint",action = "store_true",default=False,help="Do not print any detail results.")
        parser.add_argument('-l',"--list_tissue",action = "store_true",default = False,help="List tissue names in database.")
        return parser

    def run_cmd(self,args):
        args.foldchange = float(args.foldchange)
        args.weight = float(args.weight)
        args.pvalue = float(args.pvalue)
        if args.tissue.find(",") > -1:
            tissues = args.tissue
            temptis  = []
            for tis in tissues.split(","):
                temptis.append(str.capitalize(tis))
            args.tissue = ",".join(temptis)
        else:
            args.tissue = str.capitalize(args.tissue)
        args.species = str.capitalize(args.species)
        if args.species == "Mouse" and args.target == "cancersea":
            print("Error target database for mouse genome. Cancersea can't used on mouse genomes. Please use cellmarker database instead.")
            exit(0)
        if args.norefdb and not args.MarkerDB:
            print("User-defined marker database must be defined first (-M).")
            exit(0)
        rdbname = Process.check_db(args.db)
        anno = Annotator(args)
        anno.load_pickle_module(rdbname)
        outs = anno.run_detail_cmd()
        print("#Cluster","Type","Celltype","Score","Times")
        for o in outs:
            print(o)
        pass

    @staticmethod
    def list_tissue(args):
        rdbname = Process.check_db(args.db)
        anno = Annotator(args)
        anno.load_pickle_module(rdbname)
        anno.get_list_tissue("Human")
        anno.get_list_tissue("Mouse")
        sys.exit(0)

    @staticmethod
    def check_db(dbname):
        if not os.path.exists(dbname):
            dirname,fname = os.path.split(os.path.realpath(__file__))
            tempname = dirname + "/" + dbname
            if not os.path.exists(tempname):
                print(tempname)
                print("Database does not exists!",dbname)
                sys.exit(0)
        return dbname



if __name__ == "__main__":
    p = Process()
    parser = p.get_parser()
    args = parser.parse_args()
    if args.list_tissue:
        p.list_tissue(args)
    p.run_cmd(args)
