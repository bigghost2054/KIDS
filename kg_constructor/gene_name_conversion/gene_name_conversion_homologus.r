setwd("C:\\Users\\Bigghost\\Documents\\GitHub\\KIDS\\kg_constructor\\gene_name_conversion")

gene_mapping_table = read.csv("GeneMappingTable.csv",header=T,row.names=1,check.names = F)
gene_alias_table_sal = read.table("gene_result_sal.txt",header=T,check.names=F,sep="\t",quote="")
gene_alias_table_ecoli = read.table("gene_result_ecoli.txt",header=T,check.names=F,sep="\t",quote="")

target_dir = "Salm_TF"
target_file = "regulatory.txt"

input = read.table(paste(target_dir,"/",target_file,sep=""),header=T,check.names=F,sep="\t",quote="")


alias_list_sal = strsplit(gsub(" ","",gene_alias_table_sal$Aliases),",")
alias_list_sal = lapply(1:length(alias_list_sal), function(x){c(alias_list_sal[[x]],as.character(gene_alias_table_sal$Symbol[x]))})

alias_list_ecoli = strsplit(gsub(" ","",gene_alias_table_ecoli$Aliases),",")
alias_list_ecoli = lapply(1:length(alias_list_ecoli), function(x){c(alias_list_ecoli[[x]],as.character(gene_alias_table_ecoli$Symbol[x]))})

alias_list_map_sal = cbind(unlist(alias_list_sal),unlist(sapply(1:length(alias_list_sal),function(x){return(rep(x,length(alias_list_sal[[x]])))})))
alias_list_map_ecoli = cbind(unlist(alias_list_ecoli),unlist(sapply(1:length(alias_list_ecoli),function(x){return(rep(x,length(alias_list_ecoli[[x]])))})))

search_alias_idx_sal = function(gene_name){
    idx = which(alias_list_map_sal[,1] == gene_name)
    if (length(idx) > 0) return(as.numeric(alias_list_map_sal[idx[1],2]))
    return (NA)
}

search_alias_idx_ecoli = function(gene_name){
    idx = which(alias_list_map_ecoli[,1] == gene_name)
    if (length(idx) > 0) return(as.numeric(alias_list_map_ecoli[idx[1],2]))
    return (NA)
}


replace_name = function(old_names){
    new_names = old_names
    result_vec = numeric(7)
    names(result_vec) = c("InSalGeneName","InSalLocus","MapSalAliasOK","MapSalAliasFailed","MapEcoliAliasOK","MapEcoliAliasFailed","NotFound")
    
    ecoli_homologus_gene = c()
    
    for (i in 1 : length(old_names)){
        print(paste(i,"/",length(old_names)))
        #If gene in the sal table already
        if (old_names[i] %in% gene_mapping_table$SalmonellaGeneName){
            idx = which(gene_mapping_table$SalmonellaGeneName == old_names[i])
            ecoli_locus = as.character(gene_mapping_table$EcoliLocus[idx])
            alias_idx_ecoli = search_alias_idx_ecoli(ecoli_locus)
            ecoli_formal_name = as.character(gene_alias_table_ecoli$Symbol[alias_idx_ecoli])
            
            ecoli_homologus_gene = c(ecoli_homologus_gene, ecoli_formal_name)
            new_names[i] = paste(ecoli_formal_name,"(Salmonella)",sep="")
            result_vec["InSalGeneName"] = result_vec["InSalGeneName"] + 1
            next
        }
        #If gene in the sal table already
        if (old_names[i] %in% gene_mapping_table$SalmonellaLocus){
            idx = which(gene_mapping_table$SalmonellaLocus == old_names[i])
            ecoli_locus = as.character(gene_mapping_table$EcoliLocus[idx])
            alias_idx_ecoli = search_alias_idx_ecoli(ecoli_locus)
            ecoli_formal_name = as.character(gene_alias_table_ecoli$Symbol[alias_idx_ecoli])
            
            ecoli_homologus_gene = c(ecoli_homologus_gene, ecoli_formal_name)
            new_names[i] = paste(ecoli_formal_name,"(Salmonella)",sep="")
            result_vec["InSalLocus"] = result_vec["InSalLocus"] + 1
            next
        }
        
        #Check Sal Alias
        alias_idx_sal = search_alias_idx_sal(old_names[i])
        if (!is.na(alias_idx_sal)){
            sal_locus = strsplit(as.character(gene_alias_table_sal$Aliases[alias_idx_sal]),",")[[1]][1]
            
            if (sal_locus %in% gene_mapping_table$SalmonellaLocus){
                idx = which(gene_mapping_table$SalmonellaLocus == sal_locus)
                ecoli_locus = as.character(gene_mapping_table$EcoliLocus[idx])
                alias_idx_ecoli = search_alias_idx_ecoli(ecoli_locus)
                ecoli_formal_name = as.character(gene_alias_table_ecoli$Symbol[alias_idx_ecoli])
                
                ecoli_homologus_gene = c(ecoli_homologus_gene, ecoli_formal_name)
                new_names[i] = paste(ecoli_formal_name,"(Salmonella)",sep="")
                
                result_vec["MapSalAliasOK"] = result_vec["MapSalAliasOK"] + 1
            }else{
                new_names[i] = paste(old_names[i],"(Salmonella)",sep="") #NO CONVERSION
                result_vec["MapSalAliasFailed"] = result_vec["MapSalAliasFailed"] + 1
            }
            next
        }
        
        #Check Ecoli Table
        alias_idx_ecoli = search_alias_idx_ecoli(old_names[i])
        if (!is.na(alias_idx_ecoli)){
            ecoli_formal_name = as.character(gene_alias_table_ecoli$Symbol[alias_idx_ecoli])
            ecoli_locus_name = strsplit(as.character(gene_alias_table_ecoli$Aliases[alias_idx_ecoli]),",")[[1]][1]
            if (ecoli_locus_name %in% gene_mapping_table$EcoliLocus){
                ecoli_homologus_gene = c(ecoli_homologus_gene, ecoli_formal_name)
                new_names[i] = paste(ecoli_formal_name,"(Salmonella)",sep="")
                
                result_vec["MapEcoliAliasOK"] = result_vec["MapEcoliAliasOK"] + 1
            }else{
                new_names[i] = paste(old_names[i],"(Salmonella)",sep="") #NO CONVERSION
                result_vec["MapEcoliAliasFailed"] = result_vec["MapEcoliAliasFailed"] + 1
            }
            next
        }
        
        result_vec["NotFound"] = result_vec["NotFound"] + 1
    }
    return(list(new_names=new_names, ecoli_homologus_gene=ecoli_homologus_gene))
}


#Replace Subject
old_names = as.character(unique(input$Subject))
replace_res = replace_name(old_names)
new_names = replace_res$new_names
ecoli_homologus_gene = replace_res$ecoli_homologus_gene
write.csv(cbind(old_names,new_names),paste(target_dir,"/","conversion_table_subject.csv",sep=""),row.names = F)

output = input
old_subject = as.character(input$Subject)
new_subject = as.character(input$Subject)
for (i in 1 : length(old_names)){
    idx = which(old_subject == old_names[i])
    new_subject[idx] = new_names[i]
}
output$Subject = new_subject

#Replace Object
old_names = as.character(unique(input$Object))
replace_res = replace_name(old_names)
new_names = replace_res$new_names
ecoli_homologus_gene = unique(c(ecoli_homologus_gene, replace_res$ecoli_homologus_gene))
write.csv(cbind(old_names,new_names),paste(target_dir,"/","conversion_table_object.csv",sep=""),row.names = F)

old_object = as.character(input$Object)
new_object = as.character(input$Object)
for (i in 1 : length(old_names)){
    idx = which(old_object == old_names[i])
    new_object[idx] = new_names[i]
}
output$Object = new_object
write.table(output,paste(target_dir,"/",gsub(".txt","_converted.txt",target_file),sep=""),sep="\t",row.names=F,quote=F)

ecoli_homologus_gene = cbind(ecoli_homologus_gene,"is homologous to",paste(ecoli_homologus_gene,"(Salmonella)",sep=""))
colnames(ecoli_homologus_gene) = c("Subject","Predicate","Object")
write.table(ecoli_homologus_gene, paste(target_dir,"/","homologous_gene_name.txt",sep=""),row.names=F,sep="\t",quote=F)

