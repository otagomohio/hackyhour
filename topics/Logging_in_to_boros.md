# Logging in to boros

## Set up a config file (optional)

In order to speed up log in, it is a good idea to have a config file in a hidden folder called .ssh on your home drive. 

Create a plain text file called 'config' with the following format (indent with four spaces):

    Host boros
        HostName boros.otago.ac.nz
        Port 22
        User hugh

Change the text after User to your username. 

Then check if you have a folder called .ssh on your home drive:

```
ls -a
```

if needed, create the directory:

```
mkdir .ssh
```

Then put the config file there:

```
mv config .ssh/
```

## Logging in to boros from campus computer:

With config file:

```
ssh boros
```

Without config:

```
ssh user@boros.otago.ac.nz
```

## Logging in to boros from off campus:

First you have to set up VPN access. See this webpage for help:

[**Uni Otago VPN help**](https://www.otago.ac.nz/its/services/network/otago038027.html)

Once you have VPN access set up, then you can log in as above.

## Passwordless SSH login (OSX/macOS & GNU/Linux)

#### 1. Generate a key pair

On your local machine in a command prompt, as your regular user enter the following command:

```
ssh-keygen
```

Just hit Return to all of the prompts to accept the defaults. (Alternatively you can specify a passphrase that will need to be entered each time this key is used.). If you already have a keypair, ssh-keygen will complain that the id_dsa already exists, in which case you can cancel the process, and proceed to step 2.

```
    Generating public/private rsa key pair.
    Enter file in which to save the key (/home/user/.ssh/id_rsa):
    Enter passphrase (empty for no passphrase):
    Enter same passphrase again:
    Your identification has been saved in /home/user/.ssh/id_rsa.
    Your public key has been saved in /home/user/.ssh/id_rsa.pub.
    ...
```


#### 2. Copy the generated public key to remote host

```
ssh-copy-id -i ~/.ssh/id_rsa.pub username@boros.otago.ac.nz
```

The ssh-copy-id command will install the public key into your home directory on Boros (creating ~/.ssh and~/.ssh/authorized_keys with the right permissions if necessary). That's it, you're done!

If the ssh-copy-id command is not available, you'll have to manually transfer the contents of the id_rsa.pub file to ~/.ssh/authorized_keys. Here is a way to do it from your local machine:

```
cat ~/.ssh/id_rsa.pub | ssh user@boros.otago.ac.nz "cat >> ~/.ssh/authorized_keys"
```

You should now be able to login to Boros without password.

If it isn't working yet, login (with password) on Boros and make sure the ~/.ssh folder has permissions set to 700 (drwx--------), and ~/.ssh/authorized_keys is set to 600 (-rw--------)

```
  boros:~$ chmod 600 ~/.ssh/authorized_keys
  boros:~$ chmod 700 ~/.ssh
```

