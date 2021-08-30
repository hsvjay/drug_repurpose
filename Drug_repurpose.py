#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Generates a matrix - Drugs along columns and Genes along rows.
def data_preproc():
    import os
    import pandas as pd
    import glob
    import subprocess
    subprocess.run(['/bin/bash',"shell.sh"])
    with open("list_files",'w') as f:
        for ind,fl in enumerate(sorted(glob.glob("proc*.txt"))):
            if ind<1:
                pd_result_matrix=pd.read_csv(fl,'\t',header=None)
                bool_series = pd_result_matrix[0].duplicated(keep = False)     
                pd_result_matrix=(pd_result_matrix[~bool_series])
            else:
                f_pd=pd.read_csv(fl,'\t',header=None)
                bool_series = f_pd[0].duplicated(keep = False) 
                f_pd=f_pd[~bool_series]
                pd_result_matrix=pd.merge(pd_result_matrix,f_pd[[0,1]], how='outer',on=0)
            f.write(fl+"\n")
    
    pd_result_matrix.columns=["Gene.Symbol"]+["drug_"+str(i) for i in range(1,len(glob.glob("proc*.txt"))+1)]        
    pd_result_matrix.fillna(0, inplace = True) 
    pd_result_matrix.to_csv("matrix_drugs")
    with open("final_pilot.txt","w") as f:
        for ind,col in enumerate(pd_result_matrix.columns):
            if ind>0:
                f.write('\n'.join(map(str,list(pd_result_matrix[col]))))
                f.write('\n')
    return(pd_result_matrix)


# In[2]:


def data(file_name):
    import numpy as np
    all_data = np.loadtxt(file_name, unpack=True, delimiter='\t')
    data_length = int(all_data.shape[0])
    return(all_data,data_length)

def percent_change(starting_point, current_point):
    """
    Computes the percentage difference between two points
    :return: The percentage change between starting_point and current_point
    """
    default_change = 0.00001
    try:
        change = ((float(current_point) - starting_point) / abs(starting_point)) * 100.00
        if change == 0.0:
            return default_change
        else:
            return change
    except:
        return default_change

def drug_sim_fls():
    flg=0
    counter=1
    import numpy as np
    import pandas as pd
    for out_pos in range(flg,data_length-end_point,end_point):
        print(out_pos,out_pos+end_point)
    #pattern_array=[]
        with open(str(out_pos)+"_"+str(dots_for_pattern)+".txt",'w') as f:
            f.write("Gene_Start_End"+","+",".join('Drug_'+str(counter)+'_'+str(x) for x in range(counter+1,int(len(all_data)/end_point)+1))+"\n")
            for each_drug in range(out_pos,out_pos+end_point):
                pattern=[]
                if each_drug<out_pos+end_point-dots_for_pattern:
                    for index in reversed(range(dots_for_pattern)):
                        point = percent_change(all_data[each_drug], all_data[each_drug + dots_for_pattern - index])
                        pattern.append(point)
                #pattern_array.append(pattern)
                #print("inner:",each_drug)
                #print(int(len(all_data)/end_point))
                    pattern_for_all_rec=[]
            
                    for comp_drug in range(each_drug+end_point,len(all_data),end_point):
                        y=comp_drug
                    #print(y)   
                        pattern_for_rec=[]
                        for index in reversed(range(dots_for_pattern)):
                            point_r = percent_change(all_data[y], all_data[y + dots_for_pattern - index])
                            pattern_for_rec.append(point_r)

                        pattern_for_all_rec.append(pattern_for_rec)
            #print(pattern)
            #print(pattern_for_all_rec)
            #similarity
                
                    similarities_array_all=[]
                    for pat_r in pattern_for_all_rec:
                        similarities_array=[]
                        for index in range(dots_for_pattern):
            # Compute the values of similarity only if it's the first value to be computed, or if the previous one was
            # at least 50% similar
                            similarities_array.append(100.00 - abs(percent_change(pat_r[index], pattern[index])))
                        if (np.sum(np.asarray(similarities_array)>=50)>4):
                            if np.sum(np.asarray(pat_r)!=100)>5:
                                similarities_array_all.append(np.sum(similarities_array) / dots_for_pattern)
                            else:
                                similarities_array_all.append(0)
                        else:
                            similarities_array_all.append(0)
                    f.write(str(each_drug)+"_"+str(abs(each_drug+dots_for_pattern))+","+",".join(str(sim) for sim in similarities_array_all)+"\n")
        counter+=1 

def similar_drugs():
    from scipy import stats
    dic={}
    counter=1
    import pandas as pd
    for i in range(0,data_length-end_point,end_point):
        print(str(i)+"_"+str(dots_for_pattern)+".txt")
        df=pd.read_csv(str(i)+"_"+str(dots_for_pattern)+".txt")
        df=df.rank(axis=1,ascending=False,method='dense')
        ls=[]
        for ind,val in enumerate(df.columns):
            ls.append(stats.kendalltau(np.insert((np.ones(df.shape[0])),0,0),np.insert(list(df[val]),0,0))[0])
        dic["drug_"+str(counter)]=np.concatenate((np.zeros(counter),ls))
        counter+=1
    dic["drug_"+str(counter)]=np.zeros(counter)
    heat=pd.DataFrame(dic)
    tran=heat.transpose()
    heat.index=heat.columns
    tran.columns=tran.index
    htmp=heat.add(tran)
    Each_drg_max=pd.concat([heat.idxmax(),heat.max()],axis=1)
    print(Each_drg_max.nlargest(5,1))

def pattern(d1,d2):
    ele=list(all_data[(d1-1)*end_point:d1*end_point])+list(all_data[(d2-1)*end_point:d2*end_point])
    pattern_array=[]
    recog_pattern=[]
    for i in range(0,end_point-dots_for_pattern):
        print(i,i+end_point)
        pattern=[]
        rec_pattern=[]
        for index in reversed(range(dots_for_pattern)):
            point = percent_change(ele[i], ele[i + dots_for_pattern - index])
            point_r = percent_change(ele[i+end_point], ele[i+end_point + dots_for_pattern - index])
            pattern.append(point)
            rec_pattern.append(point_r)
        pattern_array.append(pattern)
        recog_pattern.append(rec_pattern)

    similarity_array=[]
    index_similarity=[]
    dic={}
    for ln in range(len(pattern_array)):
        similarity=[]
        for index in range(dots_for_pattern):
            similarity.append(100.00 - abs(percent_change(recog_pattern[ln][index], pattern_array[ln][index])))
        similarity_array.append(similarity)
        if (np.sum(np.asarray(similarity)>=50)>5):
            if np.sum(np.asarray(recog_pattern[ln])!=100)>5:
                print(ln,similarity,pattern_array[ln],recog_pattern[ln])
                dic[ln]=[similarity,pattern_array[ln],recog_pattern[ln], np.sum(similarity) / dots_for_pattern]


                
    with open("list_files",'r') as f:
        for d, val_d in enumerate(f):
            if d+1 == int(d1):
                drug1=val_d.split("_")[1]
            elif d+1 == int(d2):
                drug2=val_d.split("_")[1]
                
    genelist=list((pd.read_csv("matrix_drugs"))['Gene.Symbol'])
    import matplotlib.pyplot as plt
    for key,val in dic.items():
        xp = np.arange(0, dots_for_pattern, 1)
        plt.figure(figsize=(10, 4))
        plt.plot(genelist[key:key+dots_for_pattern], val[1], 'red', linewidth=5)
        plt.plot(genelist[key:key+dots_for_pattern], val[2], '#24BC00', linewidth=5)
        plt.xlabel("Genes",fontsize=16)
        plt.ylabel("Score",fontsize=16)
        plt.legend([drug1,drug2], loc='upper right',prop={'size': 8})
        plt.savefig(str(val[3])+".png", dpi = 72)
        print(val[3])

    kegg=[]
    for i in dic.keys():
        kegg=kegg+genelist[i:i+10]

    with open("go.txt",'w') as f:
        f.write("\n".join(item for item in kegg))



# In[4]:


if __name__=="__main__":
    import numpy as np
    import pandas as pd
    dots_for_pattern=int(input("Enter the window length for genes\n"))
    end_point=data_preproc().shape[0]
    all_data,data_length=data("final_pilot.txt")
    drug_sim_fls()
    similar_drugs()
    drug1=int(input("Enter drug - No 1\n"))
    drug2= int(input("Enter drug - No 2 \n"))
    pattern(drug1,drug2)






