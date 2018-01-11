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
  cd /data/kaijudb

  echo "downloading: http://kaiju.binf.ku.dk/database/kaiju_index.tgz"
  mkdir -p /data/kaijudb/kaiju_index
  cd /data/kaijudb/kaiju_index
  wget -q http://kaiju.binf.ku.dk/database/kaiju_index.tgz
  tar -xzf kaiju_index.tgz
  rm kaiju_index.tgz

#  echo "downloading: http://kaiju.binf.ku.dk/database/kaiju_index_pg.tgz"
#  mkdir -p /data/kaijudb/kaiju_index_pg
#  cd /data/kaijudb/kaiju_index_pg
#  wget -q http://kaiju.binf.ku.dk/database/kaiju_index_pg.tgz
#  tar -xzf kaiju_index_pg.tgz
#  rm kaiju_index_pg.tgz

#  echo "downloading: http://kaiju.binf.ku.dk/database/kaiju_index_nr.tgz"
#  mkdir -p /data/kaijudb/kaiju_index_nr
#  cd /data/kaijudb/kaiju_index_nr
#  wget -q http://kaiju.binf.ku.dk/database/kaiju_index_nr.tgz
#  tar -xzf kaiju_index_nr.tgz
#  rm kaiju_index_nr.tgz

#  echo "downloading: http://kaiju.binf.ku.dk/database/kaiju_index_nr_euk.tgz"
#  mkdir -p /data/kaijudb/kaiju_index_nr_euk
#  cd /data/kaijudb/kaiju_index_nr_euk
#  wget -q http://kaiju.binf.ku.dk/database/kaiju_index_nr_euk.tgz
#  tar -xzf kaiju_index_nr_euk.tgz
#  rm kaiju_index_nr_euk.tgz

  cd /data/kaijudb
  if [ -s "/data/kaijudb/kaiju_index/kaiju_db.fmi" ] ; then
#  if [ -s "/data/kaijudb/kaiju_index/kaiju_db.fmi" -eq 15851683711 -a -s "/data/kaijudb/kaiju_index_pg/kaiju_db.fmi" -eq 12737853741 -a -s "/data/kaijudb/kaiju_index_nr/kaiju_db_nr.fmi" -eq 45220614796 -a -s "/data/kaijudb/kaiju_index_nr_euk/kaiju_db_nr_euk.fmi" -eq 51503967892 ] ; then
#  if [ -s "/data/kaijudb/kaiju_index/kaiju_db.fmi" -a -s "/data/kaijudb/kaiju_index_pg/kaiju_db.fmi" -a -s "/data/kaijudb/kaiju_index_nr/kaiju_db_nr.fmi" -a -s "/data/kaijudb/kaiju_index_nr_euk/kaiju_db_nr_euk.fmi" ] ; then
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
