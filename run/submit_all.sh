#!/bin/bash

# List of your job folders
folders=(
  "jobs_data_2022"
  "jobs_data_2022EE"
  # "jobs_data_2023"
  # "jobs_data_2023BPix"
  "jobs_mc_2022"
  "jobs_mc_2022EE"
  # "jobs_mc_2023"
  # "jobs_mc_2023BPix"
)

# Loop through each folder
for dir in "${folders[@]}"; do
  echo "Submitting job in: $dir"
  if [ -d "$dir" ]; then
    cd "$dir"
    if [ -f "submit.sh" ]; then
      condor_submit submit.sh
    else
      echo "submit.sh not found in $dir"
    fi
    cd ..
  else
    echo "Directory $dir does not exist"
  fi
done