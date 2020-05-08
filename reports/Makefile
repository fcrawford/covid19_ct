
default: CT_extra.pdf CT_demo.pdf CT_report.pdf


CT_report.pdf: CT_report.Rmd  
	Rscript --vanilla -e "rmarkdown::render(\"CT_report.Rmd\", output_format=\"pdf_document\", output_file=\"CT_report.pdf\")"

#  Extra parameter doc from Olya
CT_extra.pdf: CT_extra.Rmd  
	Rscript --vanilla -e "rmarkdown::render(\"CT_extra.Rmd\", output_format=\"pdf_document\", output_file=\"CT_extra.pdf\")"

#  Testing visualization doc from Richard
CT_demo.pdf: CT_demo.Rmd  
	Rscript --vanilla -e "rmarkdown::render(\"CT_demo.Rmd\", output_format=\"pdf_document\", output_file=\"CT_demo.pdf\")"		

#  Compare two scenarios, not used right now, from Richard
CT_compare_report.pdf: CT_compare_report.Rmd  
	Rscript --vanilla -e "rmarkdown::render(\"CT_compare_report.Rmd\", output_format=\"pdf_document\", output_file=\"CT_compare_report.pdf\")"				

# Not used?
CT_reopen_date_report.pdf: CT_reopen_date_report.Rmd model.R run_ct_model.R params.yml
	Rscript --vanilla -e "rmarkdown::render(\"CT_reopen_date_report.Rmd\", output_format=\"pdf_document\", output_file=\"CT_reopen_date_report.pdf\")"

# Not used?
CT_testing_report.pdf: CT_testing_report.Rmd model.R run_ct_model.R params.yml
	Rscript --vanilla -e "rmarkdown::render(\"CT_testing_report.Rmd\", output_format=\"pdf_document\", output_file=\"CT_testing_report.pdf\")"

