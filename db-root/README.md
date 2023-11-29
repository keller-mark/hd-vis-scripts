
### Create the DB

```sh

srun -p interactive --pty -t 8:00:00 -n 4 --mem 16G bash
module load postgresql

cd ~/lab
mkdir postgres_sock
pg_ctl init -D $(pwd)/hd-vis-db

createdb -h $(pwd)/postgres_sock mk596
```

Modify `$(pwd)/hd-vis-db/postgresql.conf`:

```
unix_socket_directories = '/n/data1/hms/dbmi/gehlenborg/lab/postgres_sock'
```

### Start the DB server

```sh
module load postgresql
pg_ctl -D $(pwd)/hd-vis-db -l logfile start
```

### Log in to the database:

```sh
psql -h $(pwd)/postgres_sock
```

### Set a password for the database

Do this the first time starting the server and logging in

```psql
\password
```

### Use `peewee` interactively

#### On server node:

```python
from peewee import PostgresqlDatabase
psql_db = PostgresqlDatabase('mk596', host='/home/mk596/lab/postgres_sock')
psql_db.get_tables()
```

#### On other node:

Get server node hostname:

```sh
echo $(hostname)
```

```python
import os
from peewee import PostgresqlDatabase
psql_db = PostgresqlDatabase('mk596', host='compute-e-16-233.o2.rc.hms.harvard.edu', password="some_password")
psql_db.get_tables()
```


### Stop the DB server

```sh
pg_ctl -D $(pwd)/hd-vis-db stop
```
