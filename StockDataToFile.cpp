// Agron Velovic

# include <curl/curl.h>
# include <string>


void stockDataToFile(const string &tickerName,
const string &quandl_auth_token,
const string &database="WIKI",
const string &folder="quotes"){
  string mainLink="https://www.quandl.com/api/v3/datasets/";
  mainLink+=database;
  mainLink+="/"+tickerName;
  mainLink+=".csv";
  mainLink+="?sort_order=asc&auth_token=";
  mainLink+=quandl_auth_token;
  CURL *curl;
  string quandlData;
  string fName=folder;
  fName+=tickerName;
  fName+=".txt";
  curl = curl_easy_init();
  if(curl) {
  const char* linkArrayChar=mainLink.c_str();
  curl_easy_setopt(curl, CURLOPT_URL, linkArrayChar);
  curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, writerF);
  curl_easy_setopt(curl, CURLOPT_WRITEDATA, &quandlData);
  curl_easy_perform(curl);
  curl_easy_cleanup(curl);
  printToFile(fName,quandlData);
 }
}
