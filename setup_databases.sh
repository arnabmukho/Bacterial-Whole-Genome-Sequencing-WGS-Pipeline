#!/usr/bin/env bash

set -e

DBDIR=$(pwd)/db

mkdir -p "${DBDIR}"/{kraken2,checkm,antismash}

echo "======================================="
echo "Downloading CheckM database"
echo "======================================="

docker run --rm \
-v ${DBDIR}:/databases \
arnabmukho/bacterial_wgs:1.1.0-ubuntu \
bash -c "
cd /databases/checkm &&
wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz &&
tar -xzf checkm_data_2015_01_16.tar.gz
"

docker run --rm \
-v ${DBDIR}:/databases \
arnabmukho/bacterial_wgs:1.1.0-ubuntu \
checkm data setRoot /databases/checkm

echo "======================================="
echo "Downloading antiSMASH databases"
echo "======================================="

docker run --rm \
-v ${DBDIR}:/databases \
arnabmukho/bacterial_wgs:1.1.0-ubuntu \
micromamba run -n antismash_env \
download-antismash-databases \
--database-dir /databases/antismash

echo "======================================="
echo "Building Kraken2 Standard Database"
echo "======================================="

docker run --rm \
-v ${DBDIR}:/databases \
arnabmukho/bacterial_wgs:1.1.0-ubuntu \
kraken2-build \
--standard \
--threads 32 \
--db /databases/kraken2

echo "======================================="
echo "Updating Abricate Databases"
echo "======================================="

docker run --rm \
-v ${DBDIR}:/databases \
arnabmukho/bacterial_wgs:1.1.0-ubuntu \
abricate --setupdb

echo "======================================="
echo "Updating MLST Schemes"
echo "======================================="

docker run --rm \
arnabmukho/bacterial_wgs:1.1.0-ubuntu \
mlst --update

echo "======================================="
echo "Fixing Ownership"
echo "======================================="

docker run --rm \
-v ${DBDIR}:/databases \
arnabmukho/bacterial_wgs:1.1.0-ubuntu \
chown -R $(id -u):$(id -g) /databases

echo "Database installation completed."
```
