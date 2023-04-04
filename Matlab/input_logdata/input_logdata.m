cell_id = importdata('cell_id_top_plus.txt');
cell_id = string(cell_id');
gene_name = importdata('gene_name_top_plus.txt');
gene_name = string(gene_name);
zipFile = 'log_data_top_plus.txt.zip';
targetFolder = 'target';
unzip(zipFile,targetFolder);
log_data = importdata('target/log_data_top_plus.txt');
log_data = log_data.data;
save('output/mouseLiver_Slideseq_RCTD_top2000_plus.mat','-v7.3','cell_id','gene_name','log_data');


