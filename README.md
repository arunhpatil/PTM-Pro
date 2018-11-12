# PTM-Pro
<strong><h1>Post translational modifications - profiling</h1></strong>

This repository is a Command line interface bash script for Linux and Mac. This program will fetch high confident PTM sites at peptide and protein level based on user-defined ptmRS probability cut-off. A summary of the number of PTM sites that qualify the said cut-off will be generated. Also, a PTM window consisting of 7 amino acids upstream and downstream of the PTM site will be provided. This window can be used for motif enrichment or to identify conserved motifs. The PTM-window can also be used for the comparision with PhosphoSitePlus database (Nov 12, 2018). 

<strong>HOW TO USE:</strong><br>
Clone or Download the repository.<br><br>
  <strong>Usage:</strong> <br>
user@userGitHub<strong>$</strong> bash PTM_Pro.sh \<probability_cut-off> \<database> \<organism> \<input.txt><br>
  Example: bash PTM_Pro.sh 75 HsRefSeq75_NXP_cont.nr.062516.fasta human Ex_human_PTM_PSMs.txt<br><br>
  <strong>Options:</strong><br>
&nbsp;&nbsp;&nbsp;&nbsp;\<probability_cut-off>: probability cut-off ( 0 to 100)<br>
&nbsp;&nbsp;&nbsp;&nbsp;\<database>: Formatted database (Click <a href="http://ptm-pro.inhouseprotocols.com/Databases/" target="_blank">here</a> for available databases or refer the header format)<br>
&nbsp;&nbsp;&nbsp;&nbsp;\<organism>: For comparission with PhosphoSitePlus <br>
&nbsp;&nbsp;&nbsp;&nbsp;\<input.txt>: Your input text file.<br><br>
For detailed instructions for formatting and processing your input, please follow the published web portal <em><a href="http://ptm-pro.inhouseprotocols.com/" target="_blank">PTM-Pro</a></em>. <br><br>

<strong>Please Cite:</strong> Patil, A. H., Datta, K. K., et al., Dissecting Candida pathobiology: Post-translational modifications on the <em>Candida tropicalis</em> proteome. 2018. OMICS: A Journal of Integrative Biology. <a href="https://www.ncbi.nlm.nih.gov/pubmed/30106353">[PubMed]</a>
