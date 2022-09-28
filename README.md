# Standard connection process:

* If not connected from the local network, you first need to access the gate:

```$ ssh jamarino@gate.pd.infn.it```

and enter the password when prompted.
* Once in the gate or from the local network, you can connect to the t2 network:

```$ ssh t2-ui-12.pd.infn.it```

and enter your t2 password (which may differ from the gate one) when prompted.

* Enter the cms environment:

```$ cmsenv```

* and enter the proxy for gaining access to the grid:

```$ X509_USER_PROXY=/homeui/zucchett/tmp/x509up```


# About passwords:

* To change the password for jamarino@pd.infn.it (which is the same one as the gate):

https://webmail.pd.infn.it/cgi-bin/passwd.pl

* To change the password for **your** t2-network access, you need to enter, from inside the t2-ui-12, the command:

```$ passwd```

and (if I remember correctly) prompt the current password and the new password (the latter one twice).

# Copying files from remote to local:
* From the local machine:

```$ scp -r -J jamarino@gate.pd.infn.it jamarino@t2-ui-12:/homeui/jamarino/path/to/file /Users/javi/Documents/Padova/Research.nosync/ResearchActivities/```

using the option ```-r``` (recursive) to copying entire folders. The option ```-J``` first makes an scp connection to the jump host (for more info refer to the [docs](https://man7.org/linux/man-pages/man1/scp.1.html)).

* Bear in mind that the analyzed files are actually in the directory ```/lustre/cmswork/jamarino/...```. For example:

```$ scp -r -J jamarino@gate.pd.infn.it jamarino@t2-ui-12:/lustre/cmswork/jamarino/CMSSW_10_6_20/src/NanoSkim/Skimmer/OutputTest/ /Users/javi/Documents/Padova/Research.nosync/ResearchActivities/data.nosync/```




# Copying files from local to remote: (**this is still to be checked**)
* From the local machine:

```$ scp -p 2222 /path/to/the/file -J jamarino@gate.pd.infn.it jamarino@t2-ui-12:/homeui```

or 

```$ scp -P 2222 file.ext username@domain:~/```

* From remote:

```username@domain $ scp dragonmnl@local:/path/to/file.ext ~/```




