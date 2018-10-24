from flask import Flask
import requests
import re
import pandas as pd
import unicodedata
import time

app = Flask(__name__)

@app.route('/')
def index():
	return "Hello World!"
