
TTreeReaderValue<float> rv_ML_MassHH_HIGH(reader, "ML_MassHH_HIGH");
float ML_MassHH_HIGH;

int simone_vars = 1;
if (simone_vars == 1){
    ML_MassHH_HIGH = *rv_ML_MassHH_HIGH;
}