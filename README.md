# genetics
A web application about DNA sequence analysis written in Python with Flask framework.

### Functionalities
A Web-App made to analyse raw DNA sequence or read FASTA files and give various informations about its composition, pattern, melting temperature and translate it into protein sequence allowing to spot ORF easily

### How to use

#### First: Install the dependencies
This application uses flask, flask-wtf and wtforms. In order to install them you should use PIP as folowing:
```python
pip install flask, flask-wtf, wtforms
```

#### Then: modify the code depending of the way you want to deploy it
If you're running on your own computer you need to uncomment line 66 in app.py file:
```python
app.run(debug=True, host="127.0.0.1", port=5010)
```

If you're running on your a dedicated hosting solution you need to uncomment line 64 in app.py file:
```python
app.run(debug=False, host="0.0.0.0", port=5010)
```
(You can create a sub-domain that redirect you to the application via a proxy-pass using NGINX)
```nginx
    location / {
          proxy_pass http://your.ip.X.X:5010;
            }
```

To launch the application, it is very simple:
```bash
python app.py
```
(python3)



### Screenshot (with test.fasta)

#### Homepage
![Screenshot of the Application](https://i.imgur.com/hOzxmfJ.png)


#### Result Page
![Screenshot of the Application2](https://i.imgur.com/uzgX0Ro.png)
