
Create the DB:

```sh

srun -p interactive --pty -t 8:00:00 -n 4 --mem 16G bash
module load postgresql

cd ~/lab
mkdir postgres_sock
pg_ctl init -D $(pwd)/hd-vis-db

/n/app/postgresql/15.2/bin/pg_ctl -D /home/mk596/lab/hd-vis-db -l logfile start

cd ~/lab/hd-vis-db
postgres
```

Modify `$(pwd)/hd-vis-db/postgresql.conf`:

```
unix_socket_directories = '/n/data1/hms/dbmi/gehlenborg/lab/postgres_sock'
```

Start the DB server:

```sh
module load postgresql
pg_ctl -D /home/mk596/lab/hd-vis-db -l logfile start
```

Create a database:

```sh
createdb -h $(pwd)/postgres_sock mk596
```

Log in to the database:

```sh
psql -h $(pwd)/postgres_sock
```

Use `peewee` interactively:

```python
from peewee import PostgresqlDatabase
psql_db = PostgresqlDatabase('mk596', host='/home/mk596/lab/postgres_sock', user='mk596')
psql_db.get_tables()
```


Stop the DB server:

```sh
pg_ctl -D /home/mk596/lab/hd-vis-db stop
```
