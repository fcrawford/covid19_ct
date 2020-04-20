

default: CT_testing_report.pdf CT_Scenario_Report.pdf

CT_Report.pdf: CT_Report.Rmd model.R run_ct_model.R params.yml
	Rscript --vanilla -e "rmarkdown::render(\"CT_Report.Rmd\", output_format=\"pdf_document\", output_file=\"CT_Report.pdf\")"


CT_Scenario_Report.pdf: CT_Scenario_Report.Rmd model.R run_ct_model.R params.yml
	Rscript --vanilla -e "rmarkdown::render(\"CT_Scenario_Report.Rmd\", output_format=\"pdf_document\", output_file=\"CT_Scenario_Report.pdf\")"

CT_testing_report.pdf: CT_testing_report.Rmd model.R run_ct_model.R params.yml
	Rscript --vanilla -e "rmarkdown::render(\"CT_testing_report.Rmd\", output_format=\"pdf_document\", output_file=\"CT_testing_report.pdf\")"


