#!/bin/bash

. /kb/deployment/user-env.sh

python ./scripts/prepare_deploy_cfg.py ./deploy.cfg ./work/config.properties

if [ -f ./work/token ] ; then
  export KB_AUTH_TOKEN=$(<./work/token)
fi

if [ $# -eq 0 ] ; then
  sh ./scripts/start_server.sh
elif [ "${1}" = "test" ] ; then
  echo "Run Tests"
  make test
elif [ "${1}" = "async" ] ; then
  #sh ./scripts/run_async.sh
  xvfb-run bash ./scripts/run_async.sh
elif [ "${1}" = "init" ] ; then
  echo "Initialize module"
  mkdir -p /data/kaijudb
  cd /data/kaijudb

  echo "downloading: http://kaiju.binf.ku.dk/database/kaiju_db_refseq_2021-02-26.tgz"
  mkdir -p /data/kaijudb/refseq
  cd /data/kaijudb/
  wget -q http://kaiju.binf.ku.dk/database/kaiju_db_refseq_2021-02-26.tgz
  tar -xzf kaiju_db_refseq_2021-02-26.tgz
  mv kaiju_db_refseq.fmi refseq/kaiju_db_refseq.fmi
  mv names.dmp refseq/names.dmp
  mv nodes.dmp refseq/nodes.dmp
  rm kaiju_db_refseq_2021-02-26.tgz

  echo "downloading: http://kaiju.binf.ku.dk/database/kaiju_db_progenomes_2021-03-02.tgz"
  mkdir -p /data/kaijudb/progenomes
  cd /data/kaijudb/
  wget -q http://kaiju.binf.ku.dk/database/kaiju_db_progenomes_2021-03-02.tgz
  tar -xzf kaiju_db_progenomes_2021-03-02.tgz
  mv kaiju_db_progenomes.fmi progenomes/kaiju_db_progenomes.fmi
  mv names.dmp progenomes/names.dmp
  mv nodes.dmp progenomes/nodes.dmp
  rm kaiju_db_progenomes_2021-03-02.tgz

  echo "downloading: http://kaiju.binf.ku.dk/database/kaiju_db_nr_2021-02-24.tgz"
  mkdir -p /data/kaijudb/nr
  cd /data/kaijudb/
  wget -q http://kaiju.binf.ku.dk/database/kaiju_db_nr_2021-02-24.tgz
  tar -xzf kaiju_db_nr_2021-02-24.tgz
  mv kaiju_db_nr.fmi nr/kaiju_db_nr.fmi
  mv names.dmp nr/names.dmp
  mv nodes.dmp nr/nodes.dmp
  rm kaiju_db_nr_2021-02-24.tgz

  echo "downloading: http://kaiju.binf.ku.dk/database/kaiju_db_nr_euk_2021-02-24.tgz"
  mkdir -p /data/kaijudb/nr_euk
  cd /data/kaijudb/
  wget -q http://kaiju.binf.ku.dk/database/kaiju_db_nr_euk_2021-02-24.tgz
  tar -xzf kaiju_db_nr_euk_2021-02-24.tgz
  mv kaiju_db_nr_euk.fmi nr_euk/kaiju_db_nr_euk.fmi
  mv names.dmp nr_euk/names.dmp
  mv nodes.dmp nr_euk/nodes.dmp
  rm kaiju_db_nr_euk_2021-02-24.tgz

  echo "downloading: http://kaiju.binf.ku.dk/database/kaiju_db_viruses_2021-02-24.tgz"
  mkdir -p /data/kaijudb/viruses
  cd /data/kaijudb/
  wget -q http://kaiju.binf.ku.dk/database/kaiju_db_viruses_2021-02-24.tgz
  tar -xzf kaiju_db_viruses_2021-02-24.tgz
  mv kaiju_db_viruses.fmi viruses/kaiju_db_viruses.fmi
  mv names.dmp viruses/names.dmp
  mv nodes.dmp viruses/nodes.dmp
  rm kaiju_db_viruses_2021-02-24.tgz

  echo "downloading: http://kaiju.binf.ku.dk/database/kaiju_db_plasmids_2021-03-05.tgz"
  mkdir -p /data/kaijudb/plasmids
  cd /data/kaijudb/
  wget -q http://kaiju.binf.ku.dk/database/kaiju_db_plasmids_2021-03-05.tgz
  tar -xzf kaiju_db_plasmids_2021-03-05.tgz
  mv kaiju_db_plasmids.fmi plasmids/kaiju_db_plasmids.fmi
  mv names.dmp plasmids/names.dmp
  mv nodes.dmp plasmids/nodes.dmp
  rm kaiju_db_plasmids_2021-03-05.tgz

  echo "downloading: http://kaiju.binf.ku.dk/database/kaiju_db_rvdb_2021-03-05.tgz"
  mkdir -p /data/kaijudb/rvdb
  cd /data/kaijudb/
  wget -q http://kaiju.binf.ku.dk/database/kaiju_db_rvdb_2021-03-05.tgz
  tar -xzf kaiju_db_rvdb_2021-03-05.tgz
  mv kaiju_db_rvdb.fmi rvdb/kaiju_db_rvdb.fmi
  mv names.dmp rvdb/names.dmp
  mv nodes.dmp rvdb/nodes.dmp
  rm kaiju_db_rvdb_2021-03-05.tgz

  echo "downloading: http://kaiju.binf.ku.dk/database/kaiju_db_fungi_2021-03-04.tgz"
  mkdir -p /data/kaijudb/fungi
  cd /data/kaijudb/
  wget -q http://kaiju.binf.ku.dk/database/kaiju_db_fungi_2021-03-04.tgz
  tar -xzf kaiju_db_fungi_2021-03-04.tgz
  mv kaiju_db_fungi.fmi fungi/kaiju_db_fungi.fmi
  mv names.dmp fungi/names.dmp
  mv nodes.dmp fungi/nodes.dmp
  rm kaiju_db_fungi_2021-03-04.tgz

  cd /data/kaijudb

  if [ -s "/data/kaijudb/refseq/kaiju_db_refseq.fmi" -a -s "/data/kaijudb/progenomes/kaiju_db_progenomes.fmi" -a -s "/data/kaijudb/nr/kaiju_db_nr.fmi" -a -s "/data/kaijudb/nr_euk/kaiju_db_nr_euk.fmi" -a -s "/data/kaijudb/viruses/kaiju_db_viruses.fmi" -a -s "/data/kaijudb/plasmids/kaiju_db_plasmids.fmi" -a -s "/data/kaijudb/rvdb/kaiju_db_rvdb.fmi" -a -s "/data/kaijudb/fungi/kaiju_db_fungi.fmi" ] ; then
    echo "DATA DOWNLOADED SUCCESSFULLY"
    touch /data/__READY__
  else
    echo "Init failed"
  fi
elif [ "${1}" = "bash" ] ; then
  bash
elif [ "${1}" = "report" ] ; then
  export KB_SDK_COMPILE_REPORT_FILE=./work/compile_report.json
  make compile
else
  echo Unknown
fi
