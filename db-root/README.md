
```sh
srun -p interactive --pty -t 8:00:00 -n 4 --mem 16G bash
module load postgresql

cd ~/lab/hd-vis-db
postgres
```