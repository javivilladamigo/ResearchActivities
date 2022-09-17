# Standard connection process:

* If not connected from remote, you first need to access the gate:

```$ ssh jamarino@gate.pd.infn.it```

and enter the password when prompted.
* Once in the gate or from the local network, you can connect to the t2 network:

```$ ssh t2-ui-12.pd.infn.it```

and enter your t2 password (which may differ from the gate one) when prompted.

# About passwords:

* To change the password for jamarino@pd.infn.it (which is the same one as the gate):

https://webmail.pd.infn.it/cgi-bin/passwd.pl

* To change the password for **your** t2-network access, you need to enter, from inside the t2-ui-12, the command:

```$ passwd```

and (if I remember correctly) prompt the current password and the new password (the latter one twice).
