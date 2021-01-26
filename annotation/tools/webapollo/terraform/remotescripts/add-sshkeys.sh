#!/bin/bash 

# Add SSH keys

if [ -d ~/authorized-keys -a ~/.ssh ];then
    cat ~/authorized-keys/*.pub >> ~/.ssh/authorized_keys
    chmod 600 ~/.ssh/authorized_keys
    rm -rf ~/authorized-keys
fi

cat - >>$HOME/.ssh/authorized_keys <<EOF

ssh-rsa AAAAB3NzaC1yc2EAAAADAQABAAABAQCmMGgd8nvoT+uvfV+ec27D3onOkm5aFWPQH+XE6Xr62Gv12kNTgSV5l/GGbGEfkqZEc0W35LTV56WrqnFejoAeYsYbM0ZbaT7/l9PQxaNg/kC13vUVsXyYKft/dM39Nu+6exoELOS8sODBDDCxRGgSt9At+Se7lAgTnSyWFTQ++thSEzC5Ovem/146YsuuaGkTDSHv0k1Neti369fmBKhou0YeZNe3L7XSAqsgikYfUNdbd3Jv8I2j+zQ94Gk7atebewUzFqIslT1XHwSBe28g0MC68Lf3Ns0pr2mT+vdd7yiY1eQS3tvEqYuMsImL2XisvHDlTH5xm2LLucc0c2DB root@webapollo-test1

ssh-rsa AAAAB3NzaC1yc2EAAAABIwAAAgEA5YKtNXup/CiNJ0wM1rl1d+c4a4VESOND4uO1G7PYWebT+nSPjQ/nS9x1aoioiiVudBeu4gB5azQzA8Ox+JqPklaScf6b5aU6hGb5aMLp7rMGN7ACAEqGpMrnx4c4h/gx13Bs7YScD7x6+7SLTXZqrcws7yfAnd4Xg3oPFzgZ7BG4wtskyceWmi6isLZE94XMnXJeAPrClHOZcHdXyFgXDhZArrK9PLS8rQHbYxIg6fNMr49D7zPgFyqwjxG56T2UxiH73sY51bkH1OL7VdWlAwa58F2714jC5k9HQDuZMzQ66Yx5hZbB3uqq+E8K6KzfRmfZeX69eb0EjE/1hZFcZcMAddgWenyKfodvwt9iURa0PzA80xCruJMTf7Cl+5mgztFYKCuFRiMx9lfoipi0rp6QK/DnjUbXTqZaAkx+W9ufpEliGajTCXw9Q/AnSLjOAU6KRnbrVsFksGNmPC6hMwkf0Lq967H/KqZGuXhSqO057krxb8bAz73ap1xl21XOH/u9chZldL33HAEOknogJJqlbSies/e9bX3K4WkuivOMKAFAdo5ut65m9kxNv/cE7zg2cEqGZaWHZLUIvUiNi8na/HncifffB1Onhp0ABGHOCBWgaXycJx0qXYigAuDwO8wjQ/XECcVXKAbzQiGfm8gCL6dSjWCXWH9vK46QXd8= root@annotation-prod

EOF


