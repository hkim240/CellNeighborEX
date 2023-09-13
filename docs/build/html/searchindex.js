Search.setIndex({"docnames": ["api", "api/CellNeighborEX.DEanalysis.analyze_data", "api/CellNeighborEX.DEanalysis.create_nullmodel", "api/CellNeighborEX.DEanalysis.delete_files_with_keyword", "api/CellNeighborEX.DEanalysis.find_contactDEGs", "api/CellNeighborEX.DEanalysis.find_nullDEGs", "api/CellNeighborEX.DEanalysis.get_heatmap", "api/CellNeighborEX.DEanalysis.get_volcano_plot", "api/CellNeighborEX.DEanalysis.two_sample_f_test", "api/CellNeighborEX.categorization.generate_input_files", "api/CellNeighborEX.neighbors.calculate_closest_distance", "api/CellNeighborEX.neighbors.create_dataframe", "api/CellNeighborEX.neighbors.detect_neighbors", "api/CellNeighborEX.neighbors.get_neighbors", "api/CellNeighborEX.neighbors.process_dataframe", "api/CellNeighborEX.visualization.get_spatialPlot", "api/CellNeighborEX.visualization.import_expdata", "api/CellNeighborEX.visualization.set_parameters", "index", "installation", "tutorials", "tutorials/tutorial_SSV2embryo", "tutorials/tutorial_SSV2hippocampus", "tutorials/tutorial_SSV2liver_cancer", "tutorials/tutorial_seqFISH"], "filenames": ["api.rst", "api/CellNeighborEX.DEanalysis.analyze_data.rst", "api/CellNeighborEX.DEanalysis.create_nullmodel.rst", "api/CellNeighborEX.DEanalysis.delete_files_with_keyword.rst", "api/CellNeighborEX.DEanalysis.find_contactDEGs.rst", "api/CellNeighborEX.DEanalysis.find_nullDEGs.rst", "api/CellNeighborEX.DEanalysis.get_heatmap.rst", "api/CellNeighborEX.DEanalysis.get_volcano_plot.rst", "api/CellNeighborEX.DEanalysis.two_sample_f_test.rst", "api/CellNeighborEX.categorization.generate_input_files.rst", "api/CellNeighborEX.neighbors.calculate_closest_distance.rst", "api/CellNeighborEX.neighbors.create_dataframe.rst", "api/CellNeighborEX.neighbors.detect_neighbors.rst", "api/CellNeighborEX.neighbors.get_neighbors.rst", "api/CellNeighborEX.neighbors.process_dataframe.rst", "api/CellNeighborEX.visualization.get_spatialPlot.rst", "api/CellNeighborEX.visualization.import_expdata.rst", "api/CellNeighborEX.visualization.set_parameters.rst", "index.rst", "installation.rst", "tutorials.rst", "tutorials/tutorial_SSV2embryo.ipynb", "tutorials/tutorial_SSV2hippocampus.ipynb", "tutorials/tutorial_SSV2liver_cancer.ipynb", "tutorials/tutorial_seqFISH.ipynb"], "titles": ["API", "CellNeighborEX.DEanalysis.analyze_data", "CellNeighborEX.DEanalysis.create_nullmodel", "CellNeighborEX.DEanalysis.delete_files_with_keyword", "CellNeighborEX.DEanalysis.find_contactDEGs", "CellNeighborEX.DEanalysis.find_nullDEGs", "CellNeighborEX.DEanalysis.get_heatmap", "CellNeighborEX.DEanalysis.get_volcano_plot", "CellNeighborEX.DEanalysis.two_sample_f_test", "CellNeighborEX.categorization.generate_input_files", "CellNeighborEX.neighbors.calculate_closest_distance", "CellNeighborEX.neighbors.create_dataframe", "CellNeighborEX.neighbors.detect_neighbors", "CellNeighborEX.neighbors.get_neighbors", "CellNeighborEX.neighbors.process_dataframe", "CellNeighborEX.visualization.get_spatialPlot", "CellNeighborEX.visualization.import_expdata", "CellNeighborEX.visualization.set_parameters", "CellNeighborEX", "Installation", "Tutorials", "Identifying neighbor-dependent genes from Slide-seq data in a mouse embryo", "Identifying neighbor-dependent genes from Slide-seq data in mouse hippocampus", "Identifying neighbor-dependent genes from Slide-seq data in mouse liver cancer", "Identifying neighbor-dependent genes from seqFISH data in a mouse embryo"], "terms": {"import": [0, 16], "cellneighborex": [0, 19, 21, 22, 23, 24], "df_cell_id": [1, 21, 22, 23, 24], "df_gene_nam": [1, 21, 22, 23, 24], "df_log_data": [1, 21, 22, 23, 24], "path_categor": [1, 21, 22, 23, 24], "data_typ": [1, 4, 6, 7, 9, 21, 22, 23, 24], "lrcutoff": [1, 4, 5, 7, 21, 22, 23, 24], "pcutoff": [1, 4, 5, 7, 21, 22, 23, 24], "pcutoff2": [1, 4, 5, 21, 22, 23, 24], "direct": [1, 4, 5, 21, 22, 23, 24], "normality_test": [1, 4, 5, 21, 22, 23, 24], "bool": [1, 10, 14, 15], "top_gen": [1, 7, 21, 22, 23, 24], "save": [1, 6, 7, 9, 10, 14, 15, 21, 22, 23, 24], "root": [1, 9, 10, 14, 15, 21, 22, 23, 24], "de_result": [1, 21, 22, 23, 24], "function": [1, 24], "perform": [1, 4, 5], "neighbor": [1, 4, 5, 6, 18], "depend": [1, 4, 18, 19], "gene": [1, 4, 5, 6, 7, 15, 18], "express": [1, 4, 5, 15, 16, 17, 18], "analysi": [1, 9, 18], "paramet": [1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17, 21, 22, 23, 24], "datafram": [1, 7, 9, 10, 11, 14, 15, 16, 17, 21, 22, 23, 24], "contain": [1, 3, 7, 8, 9, 10, 11, 13, 14, 16, 17], "cell": [1, 2, 4, 5, 6, 7, 10, 11, 13, 14, 15, 17, 18], "id": [1, 6], "name": [1, 4, 5, 6, 7, 17, 21, 22, 23, 24], "log": [1, 2, 4, 5, 6, 7], "transform": [1, 4, 5, 24], "data": [1, 2, 4, 5, 6, 7, 8, 9, 16, 17, 18], "path": [1, 3, 16, 21, 22, 23, 24], "categor": [1, 18], "file": [1, 3, 6, 7, 9, 16], "str": [1, 15, 17, 24], "type": [1, 2, 4, 5, 6, 7, 9, 11, 12, 14, 15, 17, 18], "either": [1, 12], "imag": [1, 4, 6, 7, 9, 18, 21, 22, 23, 24], "ng": [1, 4, 6, 7, 9, 18, 21, 22, 23, 24], "float": [1, 15, 17], "ratio": [1, 4, 5, 6, 7, 21, 22, 23, 24], "cutoff": [1, 4, 5, 7, 9], "valu": [1, 4, 5, 6, 7, 8, 21, 22, 23, 24], "p": [1, 4, 5, 7, 8, 21, 22, 23, 24], "secondari": [1, 4, 5], "e": [1, 21, 22, 23, 24], "g": [1, 21, 22, 23, 24], "upregul": 1, "downregul": 1, "normal": [1, 2, 4, 5], "test": [1, 4, 5, 8, 21, 22, 23, 24], "us": [1, 4, 5, 12, 15, 18, 19, 21, 22, 23, 24], "int": 1, "number": [1, 2, 7, 9, 12, 13, 21, 22, 23, 24], "top": [1, 7, 21, 22, 23, 24], "consid": [1, 9, 12], "flag": [1, 4, 5, 10, 12, 14, 15], "indic": [1, 2, 4, 5, 6, 10, 12, 14, 15, 24], "whether": [1, 4, 5, 10, 12, 14, 15], "result": [1, 21, 22, 23, 24], "folder": [1, 6, 7, 21, 22, 23, 24], "default": [1, 9, 10, 14, 15], "i": [1, 9, 10, 14, 15, 18, 19, 21, 22, 23, 24], "return": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 21, 22, 23], "total_df": 1, "aggreg": 1, "seednumb": 2, "randsiz": 2, "prop": 2, "matchcomb": [2, 4, 5, 6], "clusterselect": [2, 4, 5, 6], "log_data": [2, 4, 5, 6], "gener": [2, 7, 9, 12, 15], "artifici": [2, 5], "heterotyp": 2, "spot": [2, 5, 9, 18, 24], "comput": [2, 18], "properti": 2, "seed": 2, "random": 2, "size": [2, 9, 15, 17, 21, 22, 23, 24], "proport": 2, "bead": [2, 17], "arrai": [2, 6, 21, 22, 23, 24], "match": [2, 5, 6], "combin": [2, 4, 5, 6], "specifi": [2, 3, 16], "cluster": [2, 4, 5, 6, 21, 22], "log_data_artificialheterospot": [2, 5], "normalized_prop": 2, "directori": [3, 9, 10, 14, 15, 16, 19, 21, 22, 23, 24], "keyword": [3, 24], "delet": [3, 24], "given": [3, 6, 9], "filenam": 3, "The": [3, 4, 5, 7, 8, 9, 14, 18, 19, 21, 22, 23, 24], "search": 3, "none": [3, 6, 7, 9, 15, 21, 22, 23, 24], "center_celltyp": [4, 5, 6, 7], "neicombuniqu": [4, 5, 6], "gene_nam": [4, 5, 6], "null_deg": [4, 5], "fdr_null": [4, 5], "logratio_nul": [4, 5], "find": [4, 5, 18], "differenti": [4, 5, 18], "deg": [4, 5, 6, 9, 21, 22, 23, 24], "relat": 4, "contact": [4, 6], "It": [4, 18], "ha": 4, "center": [4, 5, 6, 7], "comparison": 4, "list": [4, 5, 6, 10, 19], "select": [4, 5, 6, 21, 22, 23, 24], "uniqu": [4, 5, 6], "fals": [4, 5, 21, 22, 23, 24], "discoveri": [4, 5, 21, 22, 23, 24], "rate": [4, 5, 21, 22, 23, 24], "identifi": [4, 5, 18], "from": [4, 11, 16, 18, 19], "null": [4, 5], "model": 4, "fdr": [4, 5], "up": [4, 5, 18, 21, 22, 23, 24], "regul": [4, 5, 18, 21, 22, 23, 24], "down": [4, 5, 18, 21, 22, 23, 24], "filter": [4, 5], "shapiro": [4, 5], "wilk": [4, 5], "differ": 4, "tupl": [4, 8, 15, 17], "For": [4, 18, 21, 22, 23, 24], "cellcontact_degs_idx": [4, 6], "cellcontact_deg": 4, "logratio1_cellcontact": [4, 6], "1": [4, 6, 9, 21, 22, 23, 24], "pvalue1_cellcontact": 4, "fdr1_cellcontact": 4, "logratio1_tot": 4, "all": [4, 21, 22, 23, 24], "pvalue1_tot": 4, "fdr1_total": 4, "logratio2_cellcontact": [4, 6], "2": [4, 6, 15, 21, 22, 23, 24], "onli": [4, 21, 22, 23, 24], "pvalue2_cellcontact": 4, "fdr2_cellcontact": 4, "logratio_null_cellcontact": 4, "fdr_null_cellcontact": 4, "base": [5, 7, 9, 14, 17, 18, 21, 22, 23, 24], "provid": [5, 11, 15, 21, 22, 23, 24], "heterogen": [5, 9, 23], "discov": 5, "A": [5, 8, 10, 21, 22, 23], "pvalue_nul": 5, "logratio_null_tot": 5, "total": 5, "fdr_null_tot": 5, "log_data_zvalu": [6, 21, 22, 23, 24], "cell_id": 6, "foldername2": [6, 7], "plot": [6, 7, 9, 15, 21, 22, 23, 24], "heatmap": 6, "z": 6, "df": [7, 9, 10, 11, 14, 16, 21, 22, 23, 24], "volcano": [7, 21, 22, 23, 24], "input": [7, 9, 14, 21, 22, 23, 24], "choos": 7, "between": [7, 10, 24], "cneter_celltyp": 7, "annot": [7, 21, 22, 23, 24], "data1": 8, "data2": 8, "calcul": [8, 10, 13, 24], "f": 8, "statist": [8, 21, 22, 23, 24], "two": [8, 24], "sampl": [8, 9, 21, 22, 23, 24], "first": 8, "second": 8, "f_stat": 8, "p_valu": 8, "associ": 8, "sample_s": [9, 21, 22, 23, 24], "30": [9, 21, 22, 23, 24], "min_sample_s": [9, 21, 22, 23, 24], "categorized_data": [9, 21, 22, 23, 24], "paremet": 9, "minimum": 9, "requir": [9, 19, 24], "process": [9, 14, 17, 21, 22, 23, 24], "neighbor_info": [10, 14, 24], "closest": [10, 24], "distanc": [10, 24], "panda": [10, 21, 22, 23, 24], "x": [10, 21, 22, 23, 24], "y": [10, 21, 22, 23, 24], "column": [10, 11, 14, 21, 22, 23, 24], "repres": [10, 14, 24], "coordin": [10, 11, 12, 14, 24], "each": [10, 13, 14, 18], "adata": [11, 12, 24], "coord_kei": [11, 12, 24], "celltype_kei": [11, 24], "creat": [11, 19, 24], "anndata": [11, 24], "object": [11, 13, 14, 24], "dataset": [11, 12, 21, 22, 23, 24], "kei": [11, 12, 24], "access": [11, 12, 24], "inform": [11, 14, 24], "ob": [11, 24], "spatial": [11, 12, 13, 14, 15, 18], "obsm": [11, 12, 24], "construct": 11, "barcod": [11, 14, 21, 22, 23, 24], "addit": [11, 14, 24], "knn": [12, 18, 24], "radius_valu": [12, 24], "delaunai": [12, 18, 24], "detect": [12, 18, 24], "squidpi": [12, 24], "retriev": [12, 24], "connect": [12, 13, 14, 21, 22, 23, 24], "matrix": [12, 13, 14, 24], "load": [12, 24], "grid": 12, "nearest": 12, "radiu": 12, "defin": 12, "triangul": [12, 18], "spars": [12, 13, 14], "neinum": [13, 14, 24], "counter": [13, 14], "count": [13, 14, 24], "ad": [14, 24], "other": 14, "df_bg": [15, 21, 22, 23, 24], "df_red": [15, 17, 21, 22, 23, 24], "df_blue": [15, 17, 21, 22, 23, 24], "df_black": [15, 17, 21, 22, 23, 24], "label_r": [15, 21, 22, 23, 24], "label_blu": [15, 21, 22, 23, 24], "label_black": [15, 21, 22, 23, 24], "label_gen": [15, 21, 22, 23, 24], "figsiz": [15, 21, 22, 23, 24], "spatialmap": [15, 21, 22, 23, 24], "zorder_r": [15, 21, 22, 23, 24], "0": [15, 21, 22, 23, 24], "zorder_blu": [15, 21, 22, 23, 24], "zorder_black": [15, 21, 22, 23, 24], "background": [15, 17], "red": [15, 17], "blue": [15, 17], "black": [15, 17], "label": [15, 21, 22, 23, 24], "legend": 15, "colorbar": 15, "titl": 15, "draw": [15, 21, 22, 23, 24], "order": [15, 21, 22, 23, 24], "figur": [15, 18], "width": 15, "height": 15, "csv": [16, 21, 22, 23, 24], "merg": 16, "df_process": [17, 21, 22, 23, 24], "df_exp": [17, 21, 22, 23, 24], "beadsize_bg": [17, 21, 22, 23, 24], "edgecolor_bg": [17, 21, 22, 23, 24], "beadcolor_bg": [17, 21, 22, 23, 24], "beadsize_r": [17, 21, 22, 23, 24], "beadsize_blu": [17, 21, 22, 23, 24], "beadsize_black": [17, 21, 22, 23, 24], "type_r": [17, 21, 22, 23, 24], "type_blu": [17, 21, 22, 23, 24], "type_black": [17, 21, 22, 23, 24], "set": [17, 21, 22, 23, 24], "edg": 17, "color": 17, "updat": 17, "subset": 17, "github_url": [], "http": [19, 21, 22, 23, 24], "github": 19, "com": [19, 21, 22, 23], "hkim240": 19, "approach": 18, "immedi": 18, "transcriptom": 18, "st": [18, 21, 22, 23, 24], "sinlg": 18, "cellular": [18, 21, 22], "resolut": [18, 21, 22], "work": 18, "both": 18, "where": [18, 19, 21, 22, 23, 24], "exact": 18, "locat": [18, 19, 21, 22, 23, 24], "ar": [18, 19, 21, 22, 23, 24], "avail": [18, 19], "variou": 18, "algorithm": [18, 21, 22, 23], "includ": 18, "leverag": 18, "mixtur": 18, "dissect": 18, "carri": 18, "out": 18, "valid": 18, "context": 18, "below": 18, "show": 18, "workflow": 18, "python": 19, "version": [19, 21, 22, 23, 24], "3": [19, 21, 22, 23, 24], "8": [19, 21, 22, 23, 24], "11": [19, 21, 22, 23, 24], "we": [19, 21, 22, 23, 24], "recommend": 19, "conda": 19, "environ": 19, "avoid": 19, "conflict": 19, "txt": [19, 21, 22, 23, 24], "which": 19, "can": [19, 24], "download": 19, "clone": 19, "repositori": 19, "also": 19, "myenv": 19, "n": [19, 22], "10": [19, 21, 22, 23, 24], "activ": 19, "navig": 19, "Then": 19, "pip": 19, "r": 19, "pypi": 19, "pd": [21, 22, 23, 24], "o": [21, 22, 23, 24], "print": [21, 22, 23, 24], "__version__": [21, 22, 23, 24], "5": [21, 22, 23, 24], "In": [21, 22, 23, 24], "thi": [21, 22, 23, 24], "tutori": [21, 22, 23, 24], "wa": [21, 22, 23], "obtain": [21, 22, 23], "singl": [21, 22, 23], "portal": [21, 22, 23], "singlecel": [21, 22, 23], "broadinstitut": [21, 22, 23], "org": [21, 22, 23], "single_cel": [21, 22, 23], "studi": [21, 22, 23], "scp815": [21, 22], "sensit": [21, 22], "genom": [21, 22, 23], "wide": [21, 22], "profil": [21, 22], "summari": [21, 22], "pre": [21, 22, 23, 24], "format": [21, 22, 23], "ii": [21, 22, 23, 24], "regard": [21, 22, 23], "check": [21, 22, 23, 24], "your": [21, 22, 23, 24], "getcwd": [21, 22, 23, 24], "39": [21, 22, 23, 24], "user": [21, 22, 23, 24], "kimh15": [21, 22, 23, 24], "9": [21, 22, 23, 24], "make": [21, 22, 23], "exist": [21, 22, 23], "makedir": [21, 22, 23], "wget": [21, 22, 23], "figshar": [21, 22, 23], "ndownload": [21, 22, 23], "42333738": 21, "ssembryo_log_data": 21, "42333732": 21, "ssembryo_cell_id": 21, "42333735": 21, "ssembryo_gene_nam": 21, "42333702": 21, "ssembryo_rctd": 21, "42333711": 21, "ssembryo_abbrev": 21, "2023": [21, 22, 23], "09": [21, 22, 23, 24], "13": [21, 22, 23, 24], "01": [21, 22, 23, 24], "06": [21, 22, 23, 24], "resolv": [21, 22, 23], "52": [21, 22, 23, 24], "215": [21, 22, 23], "42": [21, 22, 23, 24], "80": [21, 22, 23], "50": [21, 22, 23, 24], "184": [21, 22, 23], "124": [21, 22, 23], "443": [21, 22, 23], "request": [21, 22, 23], "sent": [21, 22, 23], "await": [21, 22, 23], "respons": [21, 22, 23], "302": [21, 22, 23], "found": [21, 22, 23], "s3": [21, 22, 23], "eu": [21, 22, 23], "west": [21, 22, 23], "amazonaw": [21, 22, 23], "pfigshar": [21, 22, 23], "u": [21, 22, 23], "amz": [21, 22, 23], "aws4": [21, 22, 23], "hmac": [21, 22, 23], "sha256": [21, 22, 23], "amp": [21, 22, 23], "credenti": [21, 22, 23], "akiaiycqyoyv5jssrooa": [21, 22, 23], "20230913": [21, 22, 23], "aws4_request": [21, 22, 23], "date": [21, 22, 23], "20230913t081008z": 21, "expir": [21, 22, 23], "signedhead": [21, 22, 23], "host": [21, 22, 23], "signatur": [21, 22, 23], "66127820b6bc191d112d6a62e79d199b395906b60c151cbc319a8df9dc9fc3a5": 21, "follow": [21, 22, 23], "08": [21, 22, 23], "218": [21, 22, 23], "36": [21, 22, 23, 24], "90": [21, 22, 23], "61": [21, 22, 23], "51": [21, 22, 24], "49": [21, 22, 23], "12": [21, 22, 23, 24], "200": [21, 22, 23], "ok": [21, 22, 23], "length": [21, 22, 23, 24], "707953744": 21, "675m": 21, "text": [21, 22, 23], "plain": [21, 22, 23], "ssembryo_l": 21, "100": [21, 22, 23, 24], "gt": [21, 22, 23], "675": 21, "16m": 21, "7mb": 21, "": [21, 22, 23, 24], "67": [21, 22, 24], "17": [21, 22, 23, 24], "mb": [21, 22, 23], "20230913t081118z": 21, "75f1176fd41f4144a20b3175b3a63fd543d8dba2355f20ff828f55fdc637bf7a": 21, "18": [21, 22, 23, 24], "96": [21, 22, 23, 24], "97": [21, 22], "59": [21, 22, 24], "121": 21, "224": 21, "635430": 21, "621k": 21, "ssembryo_c": 21, "620": 21, "54k": 21, "140kb": 21, "4": [21, 22, 23, 24], "24": [21, 22, 23, 24], "126": 21, "kb": [21, 22, 23], "20230913t081125z": 21, "1245e166babff6329a4c60f62de930e0744588ad8f020e0861fab7513d3c68d3": 21, "26": [21, 22, 23, 24], "41": [21, 22, 23, 24], "211": 21, "14815": 21, "14k": [21, 23], "ssembryo_g": 21, "14": [21, 22, 23, 24], "47k": 21, "86": [21, 22], "5kb": [21, 22], "27": [21, 22, 23, 24], "20230913t081128z": 21, "027fbca56ca82c698fc23016093dc3e1bc92f3ee7349aed132ce3d7c567d04e4": 21, "28": [21, 22, 23, 24], "89": [21, 22, 23], "235": 21, "3977645": 21, "8m": 21, "ssembryo_r": 21, "79m": 21, "17mb": 21, "31": [21, 22, 23, 24], "20230913t081132z": 21, "55dd8087f41925dfe57c85cab3c6cfd2727e0c0e515867d99c028d0a78f752e3": 21, "32": [21, 22, 24], "1459": 21, "4k": 21, "ssembryo_a": 21, "42k": 21, "33": [21, 22, 23, 24], "20": [21, 22, 23, 24], "read_csv": [21, 22, 23, 24], "header": [21, 22, 23], "add": [21, 22], "abbrevi": [21, 22], "typ": [21, 22], "df_abbrev": [21, 22], "rang": [21, 22], "len": [21, 22, 23, 24], "ct1": [21, 22], "first_typ": [21, 22, 23, 24], "ct2": [21, 22], "second_typ": [21, 22, 23, 24], "j": [21, 22], "celltype1": [21, 22, 23, 24], "abbrev": [21, 22], "celltype2": [21, 22, 23, 24], "head": [21, 22, 23, 24], "prop1": [21, 22, 23, 24], "prop2": [21, 22, 23, 24], "tttttttttttttt": [21, 23], "16": [21, 22, 23, 24], "22": [21, 22, 23, 24], "sc": [21, 24], "del": 21, "5626": 21, "342930": 21, "136638": 21, "ggtttttttttttt": 21, "3738": 21, "900": 21, "3201": 21, "384132": 21, "126746": 21, "ttcaattctcgctt": 21, "my": 21, "2113": 21, "5092": 21, "7": [21, 22, 23, 24], "473911": 21, "363412": 21, "tttctgtgcagacc": 21, "2112": 21, "4926": 21, "6": [21, 22, 23, 24], "602552": 21, "000000": [21, 22, 23, 24], "tttcacagcgacgt": 21, "1789": 21, "300": 21, "4874": 21, "311685": 21, "411010": 21, "categorzi": [21, 22, 23, 24], "generate_input_fil": [21, 22, 23, 24], "15": [21, 22, 23, 24], "delimit": [21, 22, 23, 24], "t": [21, 22, 23, 24], "42352": 21, "42353": 21, "42354": 21, "42355": 21, "42356": 21, "42357": 21, "42358": 21, "42359": 21, "42360": 21, "42361": 21, "38416": 21, "78334": 21, "37038": 21, "00000": [21, 23], "42926": 21, "40604": 21, "52326": 21, "84741": 21, "82691": 21, "42685": 21, "46608": 21, "46166": 21, "35354": 21, "99487": 21, "36022": 21, "24555": 21, "00731": 21, "60155": 21, "90762": 21, "row": [21, 22, 23, 24], "42362": 21, "must": [21, 22, 23, 24], "same": [21, 22, 23, 24], "argument": [21, 22, 23, 24], "deanalysi": [21, 22, 23, 24], "analyze_data": [21, 22, 23, 24], "reguat": [21, 22, 23, 24], "true": [21, 22, 23, 24], "determin": [21, 22, 23, 24], "If": [21, 22, 23, 24], "parametr": [21, 22, 23, 24], "otherwis": [21, 22, 23, 24], "non": [21, 22, 23, 24], "when": [21, 22, 23, 24], "larger": [21, 22, 23, 24], "than": [21, 22, 23, 24], "parameter": [21, 22, 23, 24], "19": [21, 22, 23, 24], "deg_list": [21, 22, 23, 24], "72": [21, 22, 23], "00": [21, 22, 23, 24], "lt": [21, 22, 23, 24], "opt": [21, 22, 23, 24], "anaconda3": [21, 22, 23, 24], "env": [21, 22, 23, 24], "lib": [21, 22, 23, 24], "python3": [21, 22, 23, 24], "site": [21, 22, 23, 24], "py": [21, 22, 23, 24], "1131": [21, 22, 23, 24], "futurewarn": [21, 22, 23, 24], "support": [21, 22, 23, 24], "multi": [21, 22, 23, 24], "dimension": [21, 22, 23, 24], "index": [21, 22, 23, 24], "obj": [21, 22, 23, 24], "deprec": [21, 22, 23, 24], "remov": [21, 22, 23, 24], "futur": [21, 22, 23, 24], "convert": [21, 22, 23, 24], "numpi": [21, 22, 23, 24], "befor": [21, 22, 23, 24], "instead": [21, 22, 23, 24], "log_data_tot": [21, 22, 23, 24], "mean_valu": [21, 22, 23, 24], "np": [21, 22, 23, 24], "newaxi": [21, 22, 23, 24], "std_valu": [21, 22, 23, 24], "07": [21, 22, 24], "47": [21, 22, 23], "05": [21, 22, 23, 24], "40": [21, 22, 23, 24], "87": 21, "925": [21, 22, 23], "settingwithcopywarn": [21, 22, 23], "try": [21, 22, 23], "copi": [21, 22, 23], "slice": [21, 22, 23], "see": [21, 22, 23], "caveat": [21, 22, 23], "document": [21, 22, 23], "pydata": [21, 22, 23], "doc": [21, 22, 23], "stabl": [21, 22, 23, 24], "user_guid": [21, 22, 23], "html": [21, 22, 23, 24], "view": [21, 22, 23], "versu": [21, 22, 23], "largest_pvalu": [21, 22, 23], "0000000001": [21, 22, 23], "660": [21, 22, 23, 24], "userwarn": [21, 22, 23, 24], "fixedformatt": [21, 22, 23, 24], "should": [21, 22, 23, 24], "togeth": [21, 22, 23, 24], "fixedloc": [21, 22, 23, 24], "plt": [21, 22, 23, 24], "gca": [21, 22, 23, 24], "set_yticklabel": [21, 22, 23, 24], "fontsiz": [21, 22, 23, 24], "axi": [21, 22, 23, 24], "tick": [21, 22, 23, 24], "29": [21, 22, 23, 24], "56": [21, 22], "23": [21, 22, 23, 24], "93": [21, 22], "88": [21, 22, 23, 24], "58": [21, 22, 23, 24], "02": [21, 22, 23, 24], "54": [21, 22, 23, 24], "66": [21, 22, 24], "21": [21, 22, 23, 24], "55": [21, 22, 23], "37": [21, 22, 23], "43": [21, 22, 23, 24], "53": [21, 22, 23], "83": [21, 22, 23, 24], "03": [21, 22, 23, 24], "25": [21, 22, 23, 24], "91": [21, 22, 23, 24], "35": [21, 22, 23, 24], "04": [21, 22, 23, 24], "38": [21, 22, 23, 24], "45": [21, 22, 23, 24], "82": [21, 22, 24], "77": [21, 22, 23], "46": [21, 22, 23, 24], "44": [21, 22, 23, 24], "34": [21, 22, 23, 24], "48": [21, 22, 23, 24], "84": [21, 22], "57": [21, 22, 23, 24], "73": [21, 22, 23], "60": [21, 22, 24], "99": [21, 22, 23], "62": [21, 22, 23, 24], "64": [21, 22, 24], "78": [21, 22, 23], "65": [21, 22, 23], "68": [21, 22], "69": [21, 22, 23, 24], "71": [21, 22, 23, 24], "74": [21, 22], "75": [21, 22, 24], "76": [21, 22], "79": [21, 22, 24], "81": [21, 22, 23], "85": [21, 22, 23, 24], "94": [21, 22, 23], "63": [21, 22, 24], "92": [21, 22, 23, 24], "70": [21, 22, 24], "exampl": [21, 22, 23, 24], "cd24a": 21, "one": [21, 22, 23], "en": [21, 24], "l": 21, "path_select": [21, 22, 23, 24], "column_nam": [21, 22, 23, 24], "logdata": [21, 22, 23, 24], "zscore": [21, 22, 23, 24], "l_cd24a": 21, "homotypic1": [21, 22, 23, 24], "en_cd24a": 21, "homotypic2": [21, 22, 23, 24], "concat": [21, 22, 23, 24], "set_paramet": [21, 22, 23, 24], "600": [21, 22, 23], "map": [21, 22, 23, 24], "get_spatialplot": [21, 22, 23, 24], "endotheli": 21, "42333966": 22, "sshippo_log_data": 22, "42333957": 22, "sshippo_cell_id": 22, "42333960": 22, "sshippo_gene_nam": 22, "42333714": 22, "sshippo_rctd": 22, "42333708": 22, "sshippo_abbrev": 22, "20230913t084858z": 22, "e6eedaab628c8ab833bbf4a0bcbf6c3e3f8df3b094728bc52d6581992d08bf7d": 22, "176": 22, "925030656": 22, "882m": 22, "sshippo_lo": 22, "882": 22, "18m": 22, "6mb": 22, "20230913t085026z": 22, "79bcc550f701727c49f1ccb4d12f874795c2048ac4b4b03b0e18b672b8d9388a": 22, "116": 22, "208": 22, "661502": 22, "646k": 22, "sshippo_c": 22, "646": 22, "00k": 22, "142kb": 22, "122": 22, "20230913t085033z": 22, "c5df3bdb95ac577869652a4d5f934bdd8bc1664bf446ab380dd82759120aa80d": 22, "139": 22, "16604": 22, "16k": 22, "sshippo_g": 22, "21k": 22, "20230913t085035z": 22, "ee4a23be81b07c0aca09ed91055fbbd2d1c4e27ebdf1fc6cf48da654aa05cd7d": 22, "106": 22, "123": 22, "3054609": 22, "9m": 22, "sshippo_rc": 22, "91m": 22, "664kb": 22, "330": 22, "20230913t085047z": 22, "0b6197e11ae822a5c7554f82fa03bdf48769fab066b3ca2d1e9efbd3eeff1ae5": 22, "117": 22, "171": 22, "104": 22, "146": 22, "538": 22, "sshippo_ab": 22, "aacgtcataatcgt": 22, "ento": 22, "403685": 22, "103494": 22, "888": 22, "3219": 22, "tactttagcgcagt": 22, "215391": 22, "208589": 22, "4762": 22, "5020": 22, "catgcctgggttcg": 22, "412566": 22, "102875": 22, "886": 22, "3199": 22, "tcgatatggcacaa": 22, "258634": 22, "106450": 22, "2237": 22, "5144": 22, "ttatctgacgaagc": 22, "365695": 22, "132578": 22, "1031": 22, "2425": 22, "41334": 22, "41335": 22, "41336": 22, "41337": 22, "41338": 22, "41339": 22, "41340": 22, "41341": 22, "41342": 22, "41343": 22, "727628": 22, "094416": 22, "558194": 22, "295947": 22, "664255": 22, "986273": 22, "041700": 22, "330841": 22, "704539": 22, "641808": 22, "378773": 22, "259590": 22, "701731": 22, "252808": 22, "538610": 22, "729950": 22, "403333": 22, "713501": 22, "351124": 22, "110042": 22, "525093": 22, "678198": 22, "069987": 22, "41344": 22, "95": 22, "98": 22, "fabp7": 22, "ent": 22, "a_fabp7": 22, "ent_fabp7": 22, "endothelial_tip": 22, "astrocyt": 22, "scp1278": 23, "enabl": 23, "modal": 23, "clonal": 23, "tissu": 23, "42334083": 23, "ssliver_log_data": 23, "42334077": 23, "ssliver_cell_id": 23, "42334080": 23, "ssliver_gene_nam": 23, "42333705": 23, "ssliver_rctd": 23, "20230913t091214z": 23, "fe513216564b9b845c0a59c8c95f69cb7094a2892a32a24f3176c0e3594aace": 23, "144": 23, "107": 23, "412576983": 23, "393m": 23, "ssliver_lo": 23, "393": 23, "46m": 23, "2mb": 23, "20230913t091240z": 23, "b80473a58f3a5b9d112e060fb33b42a9a6fda43878149eb4c81b65c2141c3444": 23, "112": 23, "178": 23, "250": 23, "342615": 23, "335k": 23, "ssliver_c": 23, "334": 23, "58k": 23, "504kb": 23, "504": 23, "20230913t091243z": 23, "a9344167a9779e0f875b25e5b6bb38f29bddc39d5957fda820f6ab0b00ab0c06": 23, "13897": 23, "ssliver_g": 23, "57k": 23, "470": 23, "20230913t091245z": 23, "1282c94e22cab9508e0f4e016efe87a9011e759ec87d799e1e3c93f0cc67aaab": 23, "1633094": 23, "6m": 23, "ssliver_rc": 23, "56m": 23, "61mb": 23, "cttgtggttgcaga": 23, "tumoriii": 23, "monocyt": 23, "3670": 23, "2839": 23, "532488": 23, "149015": 23, "hepatocytei": 23, "2584": 23, "1702": 23, "470880": 23, "099647": 23, "acatgtctatgtta": 23, "hepatocyteii": 23, "5226": 23, "4028": 23, "737157": 23, "261249": 23, "gcgctttcgttcca": 23, "4708": 23, "1851": 23, "825412": 23, "167525": 23, "acaattgtctgcac": 23, "tumorii": 23, "1385": 23, "4176": 23, "745914": 23, "22831": 23, "22832": 23, "22833": 23, "22834": 23, "22835": 23, "22836": 23, "22837": 23, "22838": 23, "22839": 23, "22840": 23, "449814": 23, "483449": 23, "985866": 23, "943372": 23, "046288": 23, "474934": 23, "11336": 23, "250002": 23, "506298": 23, "739018": 23, "962621": 23, "22841": 23, "f13a1": 23, "monocyte_f13a1": 23, "tumoriii_f13a1": 23, "sq": 24, "scanpi": 24, "readthedoc": 24, "io": 24, "api": 24, "19416": 24, "351": 24, "create_datafram": 24, "celltype_mapped_refin": 24, "option": 24, "closest_dist": 24, "calculate_closest_dist": 24, "16it": 24, "detect_neighbor": 24, "get_neighbor": 24, "process_datafram": 24, "embryo1_pos0_cell10_z5": 24, "later": 24, "plate": 24, "mesoderm": 24, "intermedi": 24, "708437": 24, "707126": 24, "embryo1_pos0_cell100_z2": 24, "erythroid": 24, "low": 24, "qualiti": 24, "961726": 24, "943951": 24, "embryo1_pos0_cell102_z2": 24, "gut": 24, "tube": 24, "983054": 24, "284372": 24, "embryo1_pos0_cell102_z5": 24, "endothelium": 24, "993388": 24, "200491": 24, "embryo1_pos0_cell103_z5": 24, "allantoi": 24, "983855": 24, "634190": 24, "To": 24, "run": 24, "specif": 24, "need": 24, "exclud": 24, "unwant": 24, "specific_keyword": 24, "delete_files_with_keyword": 24, "neicombunique_categorized_gut": 24, "neicombunique_categorized_low": 24, "index_categorized_low": 24, "matchcomb_categorized_low": 24, "forebrain": 24, "midbrain": 24, "hindbrain": 24, "prop_categorized_surfac": 24, "ectoderm": 24, "spinal": 24, "cord": 24, "prop_categorized_low": 24, "neicombunique_categorized_forebrain": 24, "index_categorized_spin": 24, "neicombunique_categorized_splanchn": 24, "prop_categorized_forebrain": 24, "matchcomb_categorized_surfac": 24, "neicombunique_categorized_erythroid": 24, "index_categorized_surfac": 24, "index_categorized_erythroid": 24, "prop_categorized_splanchn": 24, "matchcomb_categorized_splanchn": 24, "matchcomb_categorized_erythroid": 24, "neural": 24, "crest": 24, "index_categorized_splanchn": 24, "index_categorized_forebrain": 24, "neicombunique_categorized_surfac": 24, "matchcomb_categorized_gut": 24, "prop_categorized_erythroid": 24, "neicombunique_categorized_neur": 24, "neicombunique_categorized_spin": 24, "index_categorized_gut": 24, "matchcomb_categorized_neur": 24, "prop_categorized_neur": 24, "matchcomb_categorized_forebrain": 24, "prop_categorized_spin": 24, "matchcomb_categorized_spin": 24, "prop_categorized_gut": 24, "index_categorized_neur": 24, "tolist": 24, "pp": 24, "normalize_tot": 24, "target_sum": 24, "1e4": 24, "normliz": 24, "log1p": 24, "var": 24, "to_df": 24, "reset_index": 24, "drop": 24, "embryo1_pos0_cell105_z2": 24, "embryo1_pos0_cell105_z5": 24, "embryo1_pos0_cell106_z5": 24, "embryo1_pos0_cell107_z2": 24, "embryo1_pos0_cell108_z2": 24, "embryo1_pos28_cell79_z5": 24, "embryo1_pos28_cell80_z5": 24, "embryo1_pos28_cell83_z5": 24, "embryo1_pos28_cell84_z5": 24, "embryo1_pos28_cell89_z2": 24, "embryo1_pos28_cell90_z2": 24, "embryo1_pos28_cell90_z5": 24, "embryo1_pos28_cell92_z2": 24, "embryo1_pos28_cell96_z2": 24, "embryo1_pos28_cell97_z2": 24, "241421": 24, "971863": 24, "235994": 24, "666001": 24, "313305": 24, "476862": 24, "296793": 24, "218108": 24, "660449": 24, "508516": 24, "998110": 24, "914820": 24, "003023": 24, "193126": 24, "625072": 24, "745315": 24, "865539": 24, "11797": 24, "pitx1": 24, "adjac": 24, "crest_pitx1": 24, "tube_pitx1": 24, "get": 24, "githubusercont": [], "99720939": [], "229945240": [], "2c9a2ef9": [], "2566": [], "496f": [], "9981": [], "0823cd95b813": [], "png": [], "alt": [], "align": [], "target": []}, "objects": {"CellNeighborEX": [[0, 0, 0, "-", "DEanalysis"], [0, 0, 0, "-", "categorization"], [0, 0, 0, "-", "neighbors"], [0, 0, 0, "-", "visualization"]], "CellNeighborEX.DEanalysis": [[1, 1, 1, "", "analyze_data"], [2, 1, 1, "", "create_nullmodel"], [3, 1, 1, "", "delete_files_with_keyword"], [4, 1, 1, "", "find_contactDEGs"], [5, 1, 1, "", "find_nullDEGs"], [6, 1, 1, "", "get_heatmap"], [7, 1, 1, "", "get_volcano_plot"], [8, 1, 1, "", "two_sample_f_test"]], "CellNeighborEX.categorization": [[9, 1, 1, "", "generate_input_files"]], "CellNeighborEX.neighbors": [[10, 1, 1, "", "calculate_closest_distance"], [11, 1, 1, "", "create_dataframe"], [12, 1, 1, "", "detect_neighbors"], [13, 1, 1, "", "get_neighbors"], [14, 1, 1, "", "process_dataframe"]], "CellNeighborEX.visualization": [[15, 1, 1, "", "get_spatialPlot"], [16, 1, 1, "", "import_expdata"], [17, 1, 1, "", "set_parameters"]]}, "objtypes": {"0": "py:module", "1": "py:function"}, "objnames": {"0": ["py", "module", "Python module"], "1": ["py", "function", "Python function"]}, "titleterms": {"api": 0, "immedi": 0, "neighbor": [0, 10, 11, 12, 13, 14, 21, 22, 23, 24], "detect": 0, "cell": [0, 21, 22, 23, 24], "type": [0, 21, 22, 23, 24], "categor": [0, 9, 21, 22, 23, 24], "depend": [0, 21, 22, 23, 24], "gene": [0, 21, 22, 23, 24], "express": [0, 21, 22, 23, 24], "analysi": [0, 21, 22, 23, 24], "spatial": [0, 21, 22, 23, 24], "visual": [0, 15, 16, 17, 21, 22, 23, 24], "cellneighborex": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18], "deanalysi": [1, 2, 3, 4, 5, 6, 7, 8], "analyze_data": 1, "create_nullmodel": 2, "delete_files_with_keyword": 3, "find_contactdeg": 4, "find_nulldeg": 5, "get_heatmap": 6, "get_volcano_plot": 7, "two_sample_f_test": 8, "generate_input_fil": 9, "calculate_closest_dist": 10, "create_datafram": 11, "detect_neighbor": 12, "get_neighbor": 13, "process_datafram": 14, "get_spatialplot": 15, "import_expdata": 16, "set_paramet": 17, "instal": 19, "tutori": 20, "identifi": [21, 22, 23, 24], "from": [21, 22, 23, 24], "slide": [21, 22, 23], "seq": [21, 22, 23], "data": [21, 22, 23, 24], "mous": [21, 22, 23, 24], "embryo": [21, 24], "import": [21, 22, 23, 24], "packag": [21, 22, 23, 24], "download": [21, 22, 23, 24], "load": [21, 22, 23], "bead": [21, 22, 23], "heterotyp": [21, 22, 23, 24], "spot": [21, 22, 23], "homotyp": [21, 22, 23, 24], "gener": [21, 22, 23, 24], "file": [21, 22, 23, 24], "per": [21, 22, 23, 24], "get": [21, 22, 23], "log": [21, 22, 23, 24], "normal": [21, 22, 23, 24], "perform": [21, 22, 23, 24], "context": [21, 22, 23, 24], "hippocampu": 22, "liver": 23, "cancer": 23, "seqfish": 24, "find": 24, "nearest": 24, "contact": 24, "each": 24, "produc": 24}, "envversion": {"sphinx.domains.c": 3, "sphinx.domains.changeset": 1, "sphinx.domains.citation": 1, "sphinx.domains.cpp": 9, "sphinx.domains.index": 1, "sphinx.domains.javascript": 3, "sphinx.domains.math": 2, "sphinx.domains.python": 4, "sphinx.domains.rst": 2, "sphinx.domains.std": 2, "nbsphinx": 4, "sphinx": 60}, "alltitles": {"API": [[0, "api"]], "Immediate neighbor detection": [[0, "immediate-neighbor-detection"]], "Cell type categorization": [[0, "cell-type-categorization"]], "Neighbor-dependent gene expression analysis": [[0, "neighbor-dependent-gene-expression-analysis"]], "Spatial visualization": [[0, "spatial-visualization"]], "CellNeighborEX.DEanalysis.analyze_data": [[1, "cellneighborex-deanalysis-analyze-data"]], "CellNeighborEX.DEanalysis.create_nullmodel": [[2, "cellneighborex-deanalysis-create-nullmodel"]], "CellNeighborEX.DEanalysis.delete_files_with_keyword": [[3, "cellneighborex-deanalysis-delete-files-with-keyword"]], "CellNeighborEX.DEanalysis.find_contactDEGs": [[4, "cellneighborex-deanalysis-find-contactdegs"]], "CellNeighborEX.DEanalysis.find_nullDEGs": [[5, "cellneighborex-deanalysis-find-nulldegs"]], "CellNeighborEX.DEanalysis.get_heatmap": [[6, "cellneighborex-deanalysis-get-heatmap"]], "CellNeighborEX.DEanalysis.get_volcano_plot": [[7, "cellneighborex-deanalysis-get-volcano-plot"]], "CellNeighborEX.DEanalysis.two_sample_f_test": [[8, "cellneighborex-deanalysis-two-sample-f-test"]], "CellNeighborEX.categorization.generate_input_files": [[9, "cellneighborex-categorization-generate-input-files"]], "CellNeighborEX.neighbors.calculate_closest_distance": [[10, "cellneighborex-neighbors-calculate-closest-distance"]], "CellNeighborEX.neighbors.create_dataframe": [[11, "cellneighborex-neighbors-create-dataframe"]], "CellNeighborEX.neighbors.detect_neighbors": [[12, "cellneighborex-neighbors-detect-neighbors"]], "CellNeighborEX.neighbors.get_neighbors": [[13, "cellneighborex-neighbors-get-neighbors"]], "CellNeighborEX.neighbors.process_dataframe": [[14, "cellneighborex-neighbors-process-dataframe"]], "CellNeighborEX.visualization.get_spatialPlot": [[15, "cellneighborex-visualization-get-spatialplot"]], "CellNeighborEX.visualization.import_expdata": [[16, "cellneighborex-visualization-import-expdata"]], "CellNeighborEX.visualization.set_parameters": [[17, "cellneighborex-visualization-set-parameters"]], "Installation": [[19, "installation"]], "Tutorials": [[20, "tutorials"]], "Identifying neighbor-dependent genes from Slide-seq data in a mouse embryo": [[21, "Identifying-neighbor-dependent-genes-from-Slide-seq-data-in-a-mouse-embryo"]], "Import packages": [[21, "Import-packages"], [22, "Import-packages"], [23, "Import-packages"], [24, "Import-packages"]], "Download data": [[21, "Download-data"], [22, "Download-data"], [23, "Download-data"], [24, "Download-data"]], "Load data": [[21, "Load-data"], [22, "Load-data"], [23, "Load-data"]], "Categorize Slide-seq beads into heterotypic spots and homotypic spots": [[21, "Categorize-Slide-seq-beads-into-heterotypic-spots-and-homotypic-spots"], [22, "Categorize-Slide-seq-beads-into-heterotypic-spots-and-homotypic-spots"], [23, "Categorize-Slide-seq-beads-into-heterotypic-spots-and-homotypic-spots"]], "Generate data files categorized per cell type": [[21, "Generate-data-files-categorized-per-cell-type"], [22, "Generate-data-files-categorized-per-cell-type"], [23, "Generate-data-files-categorized-per-cell-type"], [24, "Generate-data-files-categorized-per-cell-type"]], "Get log-normalized expression data": [[21, "Get-log-normalized-expression-data"], [22, "Get-log-normalized-expression-data"], [23, "Get-log-normalized-expression-data"]], "Perform neighbor-dependent gene expression analysis": [[21, "Perform-neighbor-dependent-gene-expression-analysis"], [22, "Perform-neighbor-dependent-gene-expression-analysis"], [23, "Perform-neighbor-dependent-gene-expression-analysis"], [24, "Perform-neighbor-dependent-gene-expression-analysis"]], "Visualize neighbor-dependent gene expression in the spatial context": [[21, "Visualize-neighbor-dependent-gene-expression-in-the-spatial-context"], [22, "Visualize-neighbor-dependent-gene-expression-in-the-spatial-context"], [23, "Visualize-neighbor-dependent-gene-expression-in-the-spatial-context"], [24, "Visualize-neighbor-dependent-gene-expression-in-the-spatial-context"]], "Identifying neighbor-dependent genes from Slide-seq data in mouse hippocampus": [[22, "Identifying-neighbor-dependent-genes-from-Slide-seq-data-in-mouse-hippocampus"]], "Identifying neighbor-dependent genes from Slide-seq data in mouse liver cancer": [[23, "Identifying-neighbor-dependent-genes-from-Slide-seq-data-in-mouse-liver-cancer"]], "Identifying neighbor-dependent genes from seqFISH data in a mouse embryo": [[24, "Identifying-neighbor-dependent-genes-from-seqFISH-data-in-a-mouse-embryo"]], "Find nearest neighbors for cell contact": [[24, "Find-nearest-neighbors-for-cell-contact"]], "Categorize cells into heterotypic neighbors and homotypic neighbors for each cell type": [[24, "Categorize-cells-into-heterotypic-neighbors-and-homotypic-neighbors-for-each-cell-type"]], "Produce log-normalized expression data": [[24, "Produce-log-normalized-expression-data"]], "CellNeighborEX": [[18, "cellneighborex"]]}, "indexentries": {}})