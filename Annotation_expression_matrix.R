Mean_expression = read.delim("Downloads/LVIBloodBreastCancerModel/New_cell_intensities.csv",sep=",")
Panel_file = read.delim("Downloads/LVIBloodBreastCancerModel/LVI_blood_panel.csv",sep=",")
rownames(Panel_file) = Panel_file$Metal.Tag

Order_metals = read.delim("Downloads/LVIBloodBreastCancerModel/Order_ion_channels.csv",sep=",")
Order_metals = Order_metals$X0
Order_proteins = Panel_file[Order_metals,"clean_target"]

Mean_expression = Mean_expression[,-1]
colnames(Mean_expression) = Order_proteins

pheatmap::pheatmap(cor(Mean_expression),clustering_method = 'ward')

Compensation_matrix = read.delim("Downloads/LVIBloodBreastCancerModel/compensationMatrix.csv",sep=",",row.names = 1)
Order_metals_2 = Order_metals[Order_metals%in%colnames(Compensation_matrix)]
Metals_to_add = Order_metals[!Order_metals%in%colnames(Compensation_matrix)]
Compensation_matrix = Compensation_matrix[Order_metals_2,Order_metals_2]

Compensation_matrix = cbind(Compensation_matrix,matrix(0,ncol = length(Metals_to_add),nrow =nrow(Compensation_matrix) ))
Compensation_matrix = rbind(as.matrix(Compensation_matrix),matrix(0,nrow = length(Metals_to_add),ncol =ncol(Compensation_matrix) ))
rownames(Compensation_matrix) = c(Order_metals_2,Metals_to_add)
colnames(Compensation_matrix) = c(Order_metals_2,Metals_to_add)
Compensation_matrix = Compensation_matrix[Order_metals,Order_metals]
diag(Compensation_matrix) = 1
Compensation_matrix_inverse = solve(Compensation_matrix)
Mean_expression_corrected = as.matrix(Mean_expression)%*%Compensation_matrix_inverse

colnames(Mean_expression_corrected) = Panel_file[Order_metals,"clean_target"]
pheatmap::pheatmap(cor(Mean_expression_corrected))

Compensation_matrix_2 = Compensation_matrix
diag(Compensation_matrix_2) = 0
image(Compensation_matrix_2)

write.table(Mean_expression_corrected,"Downloads/LVIBloodBreastCancerModel/Expression_matrix_annotated_corrected.csv",sep=",",row.names = FALSE)
write.table(Mean_expression,"Downloads/LVIBloodBreastCancerModel/Expression_matrix_annotated.csv",sep=",",row.names = FALSE)
