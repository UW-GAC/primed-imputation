#! /bin/bash

echo -n "-  hostname: https://imputation.biodatacatalyst.nhlbi.nih.gov
   token: " > ~/.imputationbot/imputationbot.instances

cat $1 >> ~/.imputationbot/imputationbot.instances

imputationbot instances
