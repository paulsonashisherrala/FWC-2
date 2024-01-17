#include <Arduino.h>
#ifdef ESP32
  #include <WiFi.h>
  #include <AsyncTCP.h>
#else
  #include <ESP8266WiFi.h>
  #include <ESPAsyncTCP.h>:
#endif
#include <ESPAsyncWebServer.h>
#include"matfun.h"

AsyncWebServer server(80);

const char* ssid = "lg wing";
const char* password = "password";

const char* input_parameter00 = "input00";
const char* input_parameter01 = "input01";
const char* input_parameter10 = "input10";
const char* input_parameter11 = "input11";
const char* input_parameter20 = "input20";
const char* input_parameter21 = "input21";
const char* input_parameter30 = "input30";
const char* input_parameter31 = "input31";
const char* matrix1[2]={input_parameter00,input_parameter01};     // matrix for vertex A
const char* matrix2[2]={input_parameter10,input_parameter11};     // matrix for vertex B
const char* matrix3[2]={input_parameter20,input_parameter21};     // matrix for vertex C
const char* matrix4[2]={input_parameter30,input_parameter31};     // matrix for vertex D

const char index_html[] PROGMEM = R"rawliteral(
<!DOCTYPE HTML><html><head>
    <title>TRIANGLE PROPERTIES</title>
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <style>
      html {font-family: Times New Roman; display: inline-block; text-align: center;}
      h2 {font-size: 2.0rem; color: blue;}
    </style> 
    </head><body>
    <h2>TO CHECK PROPERTIES OF TRIANGLES</h2> 
    <p>Enter the values of points A, B, C and D
    <form action="/get">
      Enter the values of Point A: <input type="number" name="input00"> <input type="number" name="input01"><br><br>
      Enter the values of Point B: <input type="number" name="input10"> <input type="number" name="input11"><br><br>
      Enter the values of Point C: <input type="number" name="input20"> <input type="number" name="input21"><br><br>
      Enter the values of Point D: <input type="number" name="input30"> <input type="number" name="input31"><br><br> 
      <input type="submit" value="Submit">
    </form><br>
  </body></html>)rawliteral";

void notFound(AsyncWebServerRequest *request) {
  request->send(404, "text/plain", "Not found");
}

void setup() {
  Serial.begin(115200);
  WiFi.mode(WIFI_STA);
  WiFi.begin(ssid, password);
  if (WiFi.waitForConnectResult() != WL_CONNECTED) {
    Serial.println("Connecting...");
    return;
  }
  Serial.println();
  Serial.print("IP Address: ");
  Serial.println(WiFi.localIP());

  server.on("/", HTTP_GET, [](AsyncWebServerRequest *request){
    request->send_P(200, "text/html", index_html);
  });

server.on("/get", HTTP_GET, [] (AsyncWebServerRequest *request) {
    double **A, **B, **C, **D, **M;
    double DM, CM, BM, AM, AC, BD, BC, CD, AB, CB;
    double CMA, BMD, DBC, BCA; // angles
    int m = 2, n = 1;
    //double l=6;length of side
    A = load_ser(request, matrix1, 2);
    B = load_ser(request, matrix2, 2);
    C = load_ser(request, matrix3, 2);
    D = load_ser(request, matrix4, 2);
    M = Matsec(A, B, m, n);

    // lengths of All sides
    DM = Matnorm(Matsub(D, M, m, n), m);
    CM = Matnorm(Matsub(C, M, m, n), m);
    BM = Matnorm(Matsub(B, M, m, n), m);
    AM = Matnorm(Matsub(A, M, m, n), m);
    AC = Matnorm(Matsub(A, C, m, n), m);
    BD = Matnorm(Matsub(B, D, m, n), m);
    BC = Matnorm(Matsub(B, C, m, n), m);
    CD = Matnorm(Matsub(C, D, m, n), m);
    AB = Matnorm(Matsub(A, B, m, n), m);
    CB = Matnorm(Matsub(C, B, m, n), m);

    // finding the angles
    CMA = angle(CM, AM, AC);
    BMD = angle(BM, DM, BD);
    DBC = angle(BC, BD, CD);
    BCA = angle(BC, AC, AB);

    String response;

    // Check conditions and add dynamic content to the response
    if ((DM == CM) && (BM == AM) && (CMA == BMD)) {
        response += "<p>&#9651; AMC &#8773; &#9651; BMD (congruent By SAS Congruency)</p>";
    } else {
        response += "<p>&#9651; AMC &#x2247; &#9651; BMD (is NOT congruent)</p>";
    }

    if (round(DBC) == 90) {
        response += "<p>&#9651; DBC is right angled At &ang; B</p>";
    } else {
        response += "<p>&#9651; DBC is NOT right angled At &ang; B</p>";
    }

    if (CB == BC && BD == AC && DBC == BCA) {
        response += "<p>&#9651; DBC &#8773; &#9651; ACB (congruent By SAS Congruency)</p>";
    } else {
        response += "<p>&#9651; DBC &#x2247; &#9651; ACB (is NOT congruent)</p>";
    }

    if (CM == AB / 2 || 2 * CM == AB) {
        response += "<p>CM = AB/2</p>";
    } else {
        response += "<p>CM &#8800; AB/2</p>";
    }
	response += "<br><a href=\"/\">Return to Home Page</a>";
    // Send the HTML response with dynamic content
    request->send(200, "text/html", response);
});
  server.onNotFound(notFound);
  server.begin();
}
void loop() { 
}
