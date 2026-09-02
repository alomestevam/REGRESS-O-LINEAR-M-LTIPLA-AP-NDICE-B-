# REGRESS-O-LINEAR-M-LTIPLA-AP-NDICE-B-
#include <WiFi.h>
#include <HTTPClient.h>
#include <Wire.h>
#include <Preferences.h>
#include <LittleFS.h>
#include <esp_sleep.h>
#include <math.h>

// ======================== REDE ========================
const char* WIFI_SSID = "SEM_SINAL";
const char* WIFI_PASSWORD = "joao3e16";
const char* SERVER_URL = "http://192.168.1.100:5000/esp32";
const uint32_t WIFI_TIMEOUT_MS = 20000UL;
const uint32_t HTTP_TIMEOUT_MS = 15000UL;

// ======================== TEMPOS ========================
const uint32_t DURACAO_EXPERIMENTO_MS = 10UL * 60UL * 1000UL;
const uint32_t INTERVALO_ENTRE_EXPERIMENTOS_S = 60UL;
const uint32_t TEMPO_ESPERA_PPK_MS = 30000UL;
const uint64_t US_POR_SEGUNDO = 1000000ULL;

// ======================== HARDWARE ========================
const uint8_t PINO_SDA = 21;
const uint8_t PINO_SCL = 22;
const uint8_t PINO_MARCADOR_EXPERIMENTO = 25;
const uint8_t PINO_MARCADOR_WIFI = 26;
const uint8_t LED_STATUS = 2;

const uint8_t MMA8452Q_ADDR = 0x1C;
const uint8_t REG_OUT_X_MSB = 0x01;
const uint8_t REG_WHO_AM_I = 0x0D;
const uint8_t REG_XYZ_DATA_CFG = 0x0E;
const uint8_t REG_CTRL_REG1 = 0x2A;

// ======================== CENÁRIOS ========================
struct ConfiguracaoCenario {
  uint8_t numero;
  uint16_t frequenciaHz;
  uint8_t intervaloMin;
};

const ConfiguracaoCenario CENARIOS[9] = {
  {1,25,1}, {2,25,3}, {3,25,5},
  {4,35,1}, {5,35,3}, {6,35,5},
  {7,45,1}, {8,45,3}, {9,45,5}
};

const uint8_t TOTAL_CENARIOS = 9;
const uint8_t REPETICOES_POR_CENARIO = 30;
const uint16_t TOTAL_EXPERIMENTOS = 270;

// ======================== MODELOS FIXOS DO PDF ========================
const double MODELO_A_INTERCEPTO = 11.32;
const double MODELO_A_COEF_F = 2.45;
const double MODELO_A_COEF_I = -0.14;

const double MODELO_B_INTERCEPTO = 25.74;
const double MODELO_B_COEF_F = 0.650;
const double MODELO_B_COEF_I = -3.96;

double preverModeloA(double F, double I) {
  return MODELO_A_INTERCEPTO + MODELO_A_COEF_F*F + MODELO_A_COEF_I*I;
}

double preverModeloB(double F, double I) {
  return MODELO_B_INTERCEPTO + MODELO_B_COEF_F*F + MODELO_B_COEF_I*I;
}

// ======================== REGRESSÃO ONLINE OLS ========================
struct AcumuladorOLS {
  uint16_t n;
  double somaF, somaI, somaFF, somaII, somaFI;
  double somaY, somaFY, somaIY, somaYY;
};

struct ResultadoOLS {
  bool valido;
  double beta0, betaF, betaI, r2, rmse;
};

AcumuladorOLS ols;

void zerarOLS() {
  memset(&ols, 0, sizeof(ols));
}

void adicionarObservacaoOLS(double F, double I, double Y) {
  ols.n++;
  ols.somaF += F; ols.somaI += I;
  ols.somaFF += F*F; ols.somaII += I*I; ols.somaFI += F*I;
  ols.somaY += Y; ols.somaFY += F*Y; ols.somaIY += I*Y;
  ols.somaYY += Y*Y;
}

bool resolver3x3(double A[3][3], double b[3], double x[3]) {
  double M[3][4];
  for (int r=0;r<3;r++) { for (int c=0;c<3;c++) M[r][c]=A[r][c]; M[r][3]=b[r]; }
  for (int p=0;p<3;p++) {
    int melhor=p;
    for (int r=p+1;r<3;r++) if (fabs(M[r][p])>fabs(M[melhor][p])) melhor=r;
    if (fabs(M[melhor][p])<1e-12) return false;
    if (melhor!=p) for (int c=p;c<4;c++) { double t=M[p][c]; M[p][c]=M[melhor][c]; M[melhor][c]=t; }
    double piv=M[p][p];
    for (int c=p;c<4;c++) M[p][c]/=piv;
    for (int r=0;r<3;r++) if (r!=p) {
      double f=M[r][p];
      for (int c=p;c<4;c++) M[r][c]-=f*M[p][c];
    }
  }
  x[0]=M[0][3]; x[1]=M[1][3]; x[2]=M[2][3];
  return true;
}

ResultadoOLS calcularOLS() {
  ResultadoOLS r{false,NAN,NAN,NAN,NAN,NAN};
  if (ols.n < 3) return r;

  double A[3][3] = {
    {(double)ols.n, ols.somaF, ols.somaI},
    {ols.somaF, ols.somaFF, ols.somaFI},
    {ols.somaI, ols.somaFI, ols.somaII}
  };
  double b[3] = {ols.somaY, ols.somaFY, ols.somaIY};
  double beta[3];
  if (!resolver3x3(A,b,beta)) return r;

  r.beta0=beta[0]; r.betaF=beta[1]; r.betaI=beta[2];
  double sse = ols.somaYY - r.beta0*ols.somaY - r.betaF*ols.somaFY - r.betaI*ols.somaIY;
  if (sse < 0 && fabs(sse) < 1e-8) sse=0;
  double sst = ols.somaYY - (ols.somaY*ols.somaY)/(double)ols.n;
  if (sst > 1e-12) r.r2 = 1.0 - sse/sst;
  r.rmse = sqrt(sse>0 ? sse/(double)ols.n : 0.0);
  r.valido=true;
  return r;
}

// ======================== PERSISTÊNCIA ========================
Preferences preferencias;
uint8_t cenarioAtual=1, repeticaoAtual=1;
bool estudoConcluido=false, estudoPausado=false;
ConfiguracaoCenario config;

void salvarOLS() { preferencias.putBytes("ols", &ols, sizeof(ols)); }
void carregarOLS() {
  if (preferencias.getBytesLength("ols") == sizeof(ols)) preferencias.getBytes("ols", &ols, sizeof(ols));
  else zerarOLS();
}

void salvarEstado() {
  preferencias.putUChar("cenario", cenarioAtual);
  preferencias.putUChar("repeticao", repeticaoAtual);
  preferencias.putBool("concluido", estudoConcluido);
  preferencias.putBool("pausado", estudoPausado);
  salvarOLS();
}

void carregarEstado() {
  preferencias.begin("exp270ols", false);
  cenarioAtual=preferencias.getUChar("cenario",1);
  repeticaoAtual=preferencias.getUChar("repeticao",1);
  estudoConcluido=preferencias.getBool("concluido",false);
  estudoPausado=preferencias.getBool("pausado",false);
  if (cenarioAtual<1 || cenarioAtual>9) cenarioAtual=1;
  if (repeticaoAtual<1 || repeticaoAtual>30) repeticaoAtual=1;
  config=CENARIOS[cenarioAtual-1];
  carregarOLS();
}

uint16_t experimentoGlobal() {
  return (uint16_t)((cenarioAtual-1)*REPETICOES_POR_CENARIO + repeticaoAtual);
}

void avancarExperimento() {
  if (repeticaoAtual < REPETICOES_POR_CENARIO) repeticaoAtual++;
  else if (cenarioAtual < TOTAL_CENARIOS) { cenarioAtual++; repeticaoAtual=1; }
  else estudoConcluido=true;
  if (!estudoConcluido) config=CENARIOS[cenarioAtual-1];
  salvarEstado();
}

void reiniciarEstudo() {
  cenarioAtual=1; repeticaoAtual=1; estudoConcluido=false; estudoPausado=false;
  config=CENARIOS[0]; zerarOLS(); salvarEstado();
}

// ======================== EXPERIMENTO ========================
struct AmostraAcelerometro { float x_g, y_g, z_g; };
uint32_t numeroAmostras=0;
double somaX=0, somaY=0, somaZ=0;
uint32_t inicioExperimentoMs=0, ultimaTransmissaoMs=0, proximaAmostraUs=0;
uint16_t transmissoesTentadas=0, transmissoesHTTP200=0;
int ultimoHTTP=0;
int32_t ultimoRSSI=0;
bool experimentoEmExecucao=false, energiaPPKRecebida=false;
double energiaPPK_mWh=NAN;

// ======================== CSV ========================
const char* CAMINHO_CSV="/PPK2_270_Experimentos_Regressao.csv";
const char* CABECALHO_CSV=
  "Cenario,Hz,Intervalo_min,Experimento,I_Medio_mA,I_Max_mA,I_Min_mA,DeltaI_mA,Energia_mWh,"
  "Previsto_Modelo_A_mWh,Previsto_Modelo_B_mWh,Residuo_Modelo_A_mWh,Residuo_Modelo_B_mWh,"
  "Amostras,Media_X_g,Media_Y_g,Media_Z_g,Transmissoes_Tentadas,Transmissoes_HTTP_200,"
  "Ultimo_HTTP_Status,Ultimo_RSSI_dBm,Duracao_ms,OLS_N,OLS_Beta0,OLS_Beta_F,OLS_Beta_I,OLS_R2,OLS_RMSE";

bool iniciarLittleFS() {
  if (!LittleFS.begin(true)) return false;
  if (!LittleFS.exists(CAMINHO_CSV)) {
    File f=LittleFS.open(CAMINHO_CSV,FILE_WRITE);
    if (!f) return false;
    f.println(CABECALHO_CSV); f.close();
  }
  return true;
}

String valorOuNA(double v, uint8_t casas=3) { return isnan(v) ? String("NA") : String(v,casas); }

void imprimirLinhaCSV(Print& out, double energiaReal, double prevA, double prevB,
                      double resA, double resB, double mediaX, double mediaY, double mediaZ,
                      uint32_t duracao, const ResultadoOLS& reg) {
  out.print("C"); out.print(config.numero); out.print(",");
  out.print(config.frequenciaHz); out.print(",");
  out.print(config.intervaloMin); out.print(",");
  out.print(repeticaoAtual); out.print(",");
  out.print("NA,NA,NA,NA,");
  out.print(valorOuNA(energiaReal)); out.print(",");
  out.print(prevA,3); out.print(","); out.print(prevB,3); out.print(",");
  out.print(valorOuNA(resA)); out.print(","); out.print(valorOuNA(resB)); out.print(",");
  out.print(numeroAmostras); out.print(",");
  out.print(mediaX,6); out.print(","); out.print(mediaY,6); out.print(","); out.print(mediaZ,6); out.print(",");
  out.print(transmissoesTentadas); out.print(","); out.print(transmissoesHTTP200); out.print(",");
  out.print(ultimoHTTP); out.print(","); out.print(ultimoRSSI); out.print(","); out.print(duracao); out.print(",");
  out.print(ols.n); out.print(",");
  out.print(reg.valido ? String(reg.beta0,6) : "NA"); out.print(",");
  out.print(reg.valido ? String(reg.betaF,6) : "NA"); out.print(",");
  out.print(reg.valido ? String(reg.betaI,6) : "NA"); out.print(",");
  out.print(reg.valido ? valorOuNA(reg.r2,6) : "NA"); out.print(",");
  out.println(reg.valido ? valorOuNA(reg.rmse,6) : "NA");
}

void gravarResultadoCSV() {
  double mediaX=numeroAmostras ? somaX/numeroAmostras : 0;
  double mediaY=numeroAmostras ? somaY/numeroAmostras : 0;
  double mediaZ=numeroAmostras ? somaZ/numeroAmostras : 0;
  double prevA=preverModeloA(config.frequenciaHz,config.intervaloMin);
  double prevB=preverModeloB(config.frequenciaHz,config.intervaloMin);
  double resA=isnan(energiaPPK_mWh)?NAN:energiaPPK_mWh-prevA;
  double resB=isnan(energiaPPK_mWh)?NAN:energiaPPK_mWh-prevB;
  ResultadoOLS reg=calcularOLS();
  uint32_t duracao=millis()-inicioExperimentoMs;

  File f=LittleFS.open(CAMINHO_CSV,FILE_APPEND);
  if (f) { imprimirLinhaCSV(f,energiaPPK_mWh,prevA,prevB,resA,resB,mediaX,mediaY,mediaZ,duracao,reg); f.close(); }

  Serial.println("===== RESULTADO_CSV =====");
  Serial.println(CABECALHO_CSV);
  imprimirLinhaCSV(Serial,energiaPPK_mWh,prevA,prevB,resA,resB,mediaX,mediaY,mediaZ,duracao,reg);
  Serial.println("===== FIM_RESULTADO_CSV =====");
}

void imprimirCSV() {
  File f=LittleFS.open(CAMINHO_CSV,FILE_READ);
  if (!f) { Serial.println("CSV_NAO_ENCONTRADO"); return; }
  Serial.println("===== INICIO_CSV_COMPLETO =====");
  while (f.available()) Serial.write(f.read());
  Serial.println("===== FIM_CSV_COMPLETO =====");
  f.close();
}

// ======================== MMA8452Q ========================
bool escreverRegistrador(uint8_t reg,uint8_t valor) {
  Wire.beginTransmission(MMA8452Q_ADDR); Wire.write(reg); Wire.write(valor);
  return Wire.endTransmission()==0;
}

bool lerRegistradores(uint8_t inicial,uint8_t* destino,size_t quantidade) {
  Wire.beginTransmission(MMA8452Q_ADDR); Wire.write(inicial);
  if (Wire.endTransmission(false)!=0) return false;
  size_t recebidos=Wire.requestFrom((int)MMA8452Q_ADDR,(int)quantidade);
  if (recebidos!=quantidade) return false;
  for (size_t i=0;i<quantidade;i++) destino[i]=Wire.read();
  return true;
}

bool iniciarMMA8452Q() {
  uint8_t id=0;
  if (!lerRegistradores(REG_WHO_AM_I,&id,1)) return false;
  Serial.printf("MMA8452Q WHO_AM_I=0x%02X\n",id);
  if (!escreverRegistrador(REG_CTRL_REG1,0x00)) return false;
  if (!escreverRegistrador(REG_XYZ_DATA_CFG,0x00)) return false;
  if (!escreverRegistrador(REG_CTRL_REG1,0x19)) return false;
  delay(10); return true;
}

int16_t converter12Bits(uint8_t msb,uint8_t lsb) {
  int16_t v=((int16_t)msb<<8)|lsb; v>>=4; if (v&0x0800) v|=0xF000; return v;
}

bool lerAcelerometro(AmostraAcelerometro& a) {
  uint8_t d[6]; if (!lerRegistradores(REG_OUT_X_MSB,d,6)) return false;
  a.x_g=converter12Bits(d[0],d[1])/1024.0f;
  a.y_g=converter12Bits(d[2],d[3])/1024.0f;
  a.z_g=converter12Bits(d[4],d[5])/1024.0f;
  return true;
}

// ======================== HTTP ========================
bool conectarWiFi() {
  WiFi.mode(WIFI_STA); WiFi.setSleep(false); WiFi.begin(WIFI_SSID,WIFI_PASSWORD);
  uint32_t ini=millis();
  while (WiFi.status()!=WL_CONNECTED) {
    if (millis()-ini>=WIFI_TIMEOUT_MS) { WiFi.disconnect(true); WiFi.mode(WIFI_OFF); ultimoHTTP=-1; return false; }
    delay(100);
  }
  ultimoRSSI=WiFi.RSSI(); return true;
}

void desligarWiFi() { WiFi.disconnect(true); WiFi.mode(WIFI_OFF); delay(50); }

String criarJSON() {
  String j; j.reserve(350);
  j+="{";
  j+="\"cenario\":"+String(config.numero)+",";
  j+="\"experimento\":"+String(repeticaoAtual)+",";
  j+="\"global\":"+String(experimentoGlobal())+",";
  j+="\"hz\":"+String(config.frequenciaHz)+",";
  j+="\"intervalo_min\":"+String(config.intervaloMin)+",";
  j+="\"amostras\":"+String(numeroAmostras)+",";
  j+="\"previsto_A_mWh\":"+String(preverModeloA(config.frequenciaHz,config.intervaloMin),3)+",";
  j+="\"previsto_B_mWh\":"+String(preverModeloB(config.frequenciaHz,config.intervaloMin),3);
  j+="}"; return j;
}

bool transmitirHTTP() {
  transmissoesTentadas++;
  digitalWrite(PINO_MARCADOR_WIFI,HIGH); digitalWrite(LED_STATUS,HIGH);
  if (!conectarWiFi()) { digitalWrite(PINO_MARCADOR_WIFI,LOW); digitalWrite(LED_STATUS,LOW); return false; }
  HTTPClient http; http.setTimeout(HTTP_TIMEOUT_MS);
  if (!http.begin(SERVER_URL)) { desligarWiFi(); digitalWrite(PINO_MARCADOR_WIFI,LOW); digitalWrite(LED_STATUS,LOW); ultimoHTTP=-2; return false; }
  http.addHeader("Content-Type","application/json");
  http.addHeader("Connection","close");
  ultimoHTTP=http.POST(criarJSON());
  if (ultimoHTTP>=200 && ultimoHTTP<300) transmissoesHTTP200++;
  Serial.printf("HTTP=%d RSSI=%ld dBm\n",ultimoHTTP,(long)ultimoRSSI);
  http.end(); desligarWiFi();
  digitalWrite(PINO_MARCADOR_WIFI,LOW); digitalWrite(LED_STATUS,LOW);
  return ultimoHTTP>=200 && ultimoHTTP<300;
}

// ======================== SERIAL / REGRESSÃO ========================
void imprimirModelos() {
  double a=preverModeloA(config.frequenciaHz,config.intervaloMin);
  double b=preverModeloB(config.frequenciaHz,config.intervaloMin);
  Serial.println("MODELO A: Y=11.32+2.45F-0.14I");
  Serial.printf("Y_A=11.32+2.45(%u)-0.14(%u)=%.3f mWh\n",config.frequenciaHz,config.intervaloMin,a);
  Serial.println("MODELO B: Y=25.74+0.650F-3.96I");
  Serial.printf("Y_B=25.74+0.650(%u)-3.96(%u)=%.3f mWh\n",config.frequenciaHz,config.intervaloMin,b);
}

void imprimirRegressaoAtual() {
  ResultadoOLS r=calcularOLS();
  Serial.println("===== REGRESSAO_MULTIPLA_ONLINE =====");
  Serial.printf("N=%u\n",ols.n);
  if (!r.valido) { Serial.println("Coeficientes ainda não identificáveis: forneça Energia_mWh para cenários com variação de F e I."); return; }
  Serial.printf("Energia_mWh = %.6f + %.6f*Hz + %.6f*Intervalo_min\n",r.beta0,r.betaF,r.betaI);
  Serial.printf("R2=%.6f RMSE=%.6f mWh\n",r.r2,r.rmse);
}

bool interpretarPPK(String comando) {
  if (!comando.startsWith("PPK,")) return false;
  String valor=comando.substring(4); valor.trim(); valor.replace(',','.');
  if (valor.equalsIgnoreCase("NA") || valor.length()==0) {
    energiaPPK_mWh=NAN; energiaPPKRecebida=true; Serial.println("PPK_ENERGIA=NA"); return true;
  }
  double y=valor.toDouble();
  if (y<=0) { Serial.println("ERRO_PPK_VALOR_INVALIDO"); return true; }
  energiaPPK_mWh=y; energiaPPKRecebida=true;
  adicionarObservacaoOLS(config.frequenciaHz,config.intervaloMin,y); salvarOLS();
  Serial.printf("PPK_ENERGIA_RECEBIDA=%.3f mWh\n",y); imprimirRegressaoAtual(); return true;
}

void processarSerial() {
  if (!Serial.available()) return;
  String c=Serial.readStringUntil('\n'); c.trim();
  if (interpretarPPK(c)) return;
  c.toUpperCase();
  if (c=="STATUS") Serial.printf("C%u R%u GLOBAL=%u/%u F=%uHz I=%umin OLS_N=%u\n",cenarioAtual,repeticaoAtual,experimentoGlobal(),TOTAL_EXPERIMENTOS,config.frequenciaHz,config.intervaloMin,ols.n);
  else if (c=="MODELOS") imprimirModelos();
  else if (c=="REGRESSAO") imprimirRegressaoAtual();
  else if (c=="CSV") imprimirCSV();
  else if (c=="PAUSAR") { estudoPausado=true; salvarEstado(); Serial.println("PAUSA_PROGRAMADA"); }
  else if (c=="RETOMAR") { estudoPausado=false; salvarEstado(); Serial.println("ESTUDO_RETOMADO"); }
  else if (c=="RESET_ESTUDO") {
    reiniciarEstudo(); if (LittleFS.exists(CAMINHO_CSV)) LittleFS.remove(CAMINHO_CSV); iniciarLittleFS(); Serial.println("ESTUDO_E_OLS_REINICIADOS");
  }
  else if (c.length()>0) Serial.println("COMANDO_DESCONHECIDO");
}

void aguardarPPK() {
  energiaPPKRecebida=false; energiaPPK_mWh=NAN;
  if (TEMPO_ESPERA_PPK_MS==0) { energiaPPKRecebida=true; return; }
  Serial.println("Informe Energia_mWh do PPK II: PPK,38.50  ou  PPK,NA");
  uint32_t ini=millis();
  while (!energiaPPKRecebida && millis()-ini<TEMPO_ESPERA_PPK_MS) { processarSerial(); delay(20); }
  if (!energiaPPKRecebida) { energiaPPK_mWh=NAN; energiaPPKRecebida=true; Serial.println("PPK_TIMEOUT -> Energia_mWh=NA"); }
}

// ======================== CICLO ========================
void pulso(uint8_t pino,uint32_t ms) { digitalWrite(pino,HIGH); delay(ms); digitalWrite(pino,LOW); }

void zerarExperimento() {
  numeroAmostras=0; somaX=somaY=somaZ=0;
  transmissoesTentadas=transmissoesHTTP200=0; ultimoHTTP=0; ultimoRSSI=0;
  energiaPPK_mWh=NAN; energiaPPKRecebida=false;
}

void entrarDeepSleep(uint32_t s) {
  digitalWrite(PINO_MARCADOR_EXPERIMENTO,LOW); digitalWrite(PINO_MARCADOR_WIFI,LOW); digitalWrite(LED_STATUS,LOW);
  desligarWiFi(); Serial.flush();
  esp_sleep_enable_timer_wakeup((uint64_t)s*US_POR_SEGUNDO);
  esp_deep_sleep_start();
}

void iniciarExperimento() {
  zerarExperimento(); experimentoEmExecucao=true;
  inicioExperimentoMs=millis(); ultimaTransmissaoMs=inicioExperimentoMs; proximaAmostraUs=micros();
  Serial.printf("INICIO,C%u,R%u,GLOBAL=%u/%u,F=%uHz,I=%umin\n",config.numero,repeticaoAtual,experimentoGlobal(),TOTAL_EXPERIMENTOS,config.frequenciaHz,config.intervaloMin);
  pulso(PINO_MARCADOR_EXPERIMENTO,300); imprimirModelos(); transmitirHTTP();
}

void finalizarExperimento() {
  experimentoEmExecucao=false;
  pulso(PINO_MARCADOR_EXPERIMENTO,150); delay(150); pulso(PINO_MARCADOR_EXPERIMENTO,150);
  Serial.printf("FIM,C%u,R%u,GLOBAL=%u,AMOSTRAS=%lu\n",config.numero,repeticaoAtual,experimentoGlobal(),(unsigned long)numeroAmostras);
  aguardarPPK(); gravarResultadoCSV(); imprimirRegressaoAtual(); avancarExperimento();
  if (estudoConcluido) {
    Serial.println("TODOS_OS_270_EXPERIMENTOS_CONCLUIDOS"); imprimirRegressaoAtual();
    while (true) { processarSerial(); delay(100); }
  }
  if (estudoPausado) { Serial.println("ESTUDO_PAUSADO"); while (estudoPausado) { processarSerial(); delay(100); } }
  entrarDeepSleep(INTERVALO_ENTRE_EXPERIMENTOS_S);
}

void setup() {
  Serial.begin(115200); Serial.setTimeout(1000); delay(800);
  pinMode(PINO_MARCADOR_EXPERIMENTO,OUTPUT); pinMode(PINO_MARCADOR_WIFI,OUTPUT); pinMode(LED_STATUS,OUTPUT);
  digitalWrite(PINO_MARCADOR_EXPERIMENTO,LOW); digitalWrite(PINO_MARCADOR_WIFI,LOW); digitalWrite(LED_STATUS,LOW);
  carregarEstado(); iniciarLittleFS(); Serial.println(CABECALHO_CSV);
  if (estudoConcluido) { Serial.println("ESTUDO_JA_CONCLUIDO"); imprimirRegressaoAtual(); while (true) { processarSerial(); delay(100); } }
  Wire.begin(PINO_SDA,PINO_SCL); Wire.setClock(400000);
  if (!iniciarMMA8452Q()) { while (true) { Serial.println("FALHA_MMA8452Q"); delay(1000); } }
  iniciarExperimento();
}

void loop() {
  processarSerial();
  if (!experimentoEmExecucao) { delay(10); return; }
  uint32_t agoraMs=millis(), agoraUs=micros();
  uint32_t periodoUs=1000000UL/config.frequenciaHz;
  uint32_t intervaloMs=(uint32_t)config.intervaloMin*60UL*1000UL;

  if ((int32_t)(agoraUs-proximaAmostraUs)>=0) {
    proximaAmostraUs=agoraUs+periodoUs;
    AmostraAcelerometro a;
    if (lerAcelerometro(a)) { somaX+=a.x_g; somaY+=a.y_g; somaZ+=a.z_g; numeroAmostras++; }
  }

  if (agoraMs-ultimaTransmissaoMs>=intervaloMs) {
    ultimaTransmissaoMs+=intervaloMs; transmitirHTTP();
  }

  if (agoraMs-inicioExperimentoMs>=DURACAO_EXPERIMENTO_MS) finalizarExperimento();
  delayMicroseconds(100);
} 
