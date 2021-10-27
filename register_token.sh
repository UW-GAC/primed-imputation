#! /bin/bash

mkdir ~/.imputationbot
printf -- "-  hostname: https://imputation.biodatacatalyst.nhlbi.nih.gov\n   token: " > ~/.imputationbot/imputationbot.instances
cat $1 >> ~/.imputationbot/imputationbot.instances

imputationbot instances
