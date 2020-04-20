

default: CT_testing_report.pdf CT_reopen_date_report.pdf



CT_reopen_date_report.pdf: CT_reopen_date_report.Rmd model.R run_ct_model.R params.yml
	Rscript --vanilla -e "rmarkdown::render(\"CT_reopen_date_report.Rmd\", output_format=\"pdf_document\", output_file=\"CT_reopen_date_report.pdf\")"

CT_testing_report.pdf: CT_testing_report.Rmd model.R run_ct_model.R params.yml
	Rscript --vanilla -e "rmarkdown::render(\"CT_testing_report.Rmd\", output_format=\"pdf_document\", output_file=\"CT_testing_report.pdf\")"


