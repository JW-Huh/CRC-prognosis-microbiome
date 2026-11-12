# 식이–장내 마이크로바이옴–대장암 예후: Letter 연구계획, 에너지 보정 방안, 매개분석 설계

작성 2026-11-12 · 허지원 (Translational Microbiome Lab, CNU) · 수행 강다영 (학부연구참여)
전제: `4-Script/run_all.R`로 2022년 논문이 수치까지 그대로 재현됨을 확인함(README 체크리스트). 이 문서는 그 위에서 무엇을, 어떤 순서로, 어떤 방법으로 할지를 정리한다.

---

## 1. 문헌 조사 요약 (2026-11-12 기준, PubMed/Europe PMC)

### 1-1. 우리 코호트(SNUH)에서 나온 논문
| 논문 | 설계 | 식이 변수 | 에너지 보정 | 예후 | 마이크로바이옴 |
|---|---|---|---|---|---|
| Huh 2022, *Microbiome* 10:203 | 333명, 수술 전 대변 16S, PFS/OS | 없음 | – | ○ (PFS/OS, M5 score) | ○ |
| Hoang T, Kim MJ, Park JW, Jeong SY, Lee J, Shin A. 2022, *BMC Cancer* 22:656 | 115명, 단면 | 106항목 SQ-FFQ → 92 영양소, 16(17) 식품군 | **residual method**(원문: "a residual method was used for the energy adjustment of the 106 food items, 16 food groups, and 92 nutrients"); 로그변환·성별 구분·0섭취 처리 언급 없음; 중앙값 이분 | ✕ (저자들이 명시적으로 "could not evaluate prognosis") | ○ (alpha/beta, LEfSe, k-medoids enterotype) |
| Hoang 2023, *BMC Microbiol* 23:33 | 331명, 단면 | 식이 다양성 트리(Procrustes) | 해당 없음 | ✕ | ○ (흡연·음주·당뇨에서 AKP 효과) |
| Jun SY … Park JW, Shin A 2024, *Cancer Med* 13:e7084 | 509명, MR | 없음(혈당 형질) | – | ○ (3y OS/PFS) | ✕ |
| Won D, … Park JW, Lee J, Cho S, Shin A 2025, *Cancer Med* 14:e71461 | 국가 검진 자료 7,553명 | **없음**(신체활동·음주·흡연·BMI) | – | ○ (사망) | ✕ |
| Han JY, Kim MJ, Park JW, Jeong SY 2026, *Ann Coloproctol* (review) | – | "dietary modulation" 소개 수준 | – | – | – |

**결론: 박지원·김민정·신애선 그룹에서 FFQ 기반 식이(영양소·식품군·식이패턴)와 CRC 예후(PFS/DFS/OS/재발)를 직접 분석한 논문은 아직 없다.** 식이–마이크로바이옴은 Hoang 2022/2023 단면 연구 두 편뿐이고, 두 편 모두 "예후는 볼 수 없었다"고 한계를 적었다. 따라서 **식이 → 수술 전 대변 마이크로바이옴 → 예후**를 잇는 분석은 그룹 내에서도, 문헌 전체에서도 공백이다.

### 1-2. 국내 비교 대상
- Kim M … Bae S (계명대) 2025, *Nutrients* 17:1522: n=126, 수술 전 109항목 FFQ → DII. **에너지 보정 없음**, DII 절단값을 OS log-rank 최대화로 선택(과적합 위험), high-DII n=28. 5y OS 41 vs 90%. 마이크로바이옴 없음. → 우리 letter에서 "에너지 보정 + 마이크로바이옴"으로 차별화할 수 있는 대조군 논문.

### 1-3. 방법론 기준이 되는 외국 코호트 (식이–CRC 생존)
| 연구 | 에너지 보정 방식(확인된 원문) | 범주화 |
|---|---|---|
| Meyerhardt 2007 *JAMA* (CALGB 89803, 1,009명) | "All nutrient values were energy-adjusted using the residuals method" + **총열량을 공변량으로 추가** | 식이패턴 5분위 |
| Song 2018 *JAMA Oncol* (NHS/HPFS, 1,575명) | "adjusted fiber intake for total caloric intake using the nutrient residual method" | 5 g/일당, 4분위 |
| Van Blarigan 2020 *JAMA Netw Open* (CALGB/SWOG 80405, 1,284명) | "adjusted for energy using the residual method" + 모델1에 kcal 포함; **성별 5분위** | 성별 5분위 |
| Van Blarigan 2018 *JAMA Oncol* (ACS score) | 점수는 절단값 기반, **총열량을 공변량**(standard multivariate) | 점수 범주 |
| van Zutphen 2021 *AJCN* (COLON/EnCoRe, 1,425명) | 절대 절단값(섬유 ≥30 g 등); 총열량은 교란요인 검정 후 제외(<5% 변화) | 1 SD당 |
| Hang 2024 *eClinicalMedicine* (UPF) | "adjusted for total energy using the residual method"(식품 servings에 residual 적용) | 5분위 |
| Cheng 2022 *JCO* (CALGB 89803 nomogram) | 에너지 보정 없음(빈도 범주) | 빈도 범주 |
| Tan 2026 *Food Funct*, Wu 2026 *JNCI* (광둥 코호트 2,799명) | isocaloric substitution + 염증지표 **매개분석**(5.6–9.4% 매개) | 4분위 |

→ 주류(CALGB/NHS 계열)의 관행: **Willett residual method를 영양소(및 servings)에 적용하고, 총열량을 Cox 공변량으로 한 번 더 넣으며, 5분위(또는 성별 5분위)로 범주화**. Hoang 2022는 residual까지는 같지만 중앙값 이분·로그변환 없음·에너지 공변량 없음.

### 1-4. 마이크로바이옴을 매개자로 본 연구
- Shimomura 2023 *Int J Cancer* 152:1752 (일본 NCC, 180 CRC vs 242 대조): 섬유 → *F. nucleatum* → CRC(발생)에 **counterfactual causal mediation**(NIE/NDE, risk-difference scale) 적용; *F. nucleatum*이 23%, *G. morbillorum*이 10.5% 매개. → **동아시아 CRC + FFQ + 종 수준 매개분석의 유일한 선례**(단, 결과는 발생이지 예후가 아님).
- Byrd 2025 *Int J Cancer* 157:64 (ColoCare, 166명, 수술 전 대변, DFS): Huh 2022의 "마이크로바이옴 → 예후" 축을 외부 재현(*Peptostreptococcus* HR 1.62, *Coprococcus* 0.49). **식이 변수 없음.**
- Deng 2026 *Br J Cancer*: pks+ *E. coli*, *F. nucleatum*이 식습관 → 조기발병 CRC를 매개(발생).
- Kim B … Lee S 2024 *Front Microbiol* 15:1282932: 한국 CRC 수술 후 환자 40명 교차설계 RCT, mMAC(microbiota-accessible carbohydrate) 식이 3주 → *Prevotella*·acetate·propionate 증가. → "식이로 Prevotella를 올릴 수 있다"는 국내 기계적 근거.
- 방법론: Luo 2020 *PLoS Comput Biol*(HIMA-survival), Zhang 2021 *Bioinformatics*(Cox 고차원 매개), Sohn & Li 2019 *Ann Appl Stat*(compositional mediation), R `HIMA`(hima_survival, compositional), `CMAverse`, `mediation`.

**공백(명시적으로 미발표):** FFQ 식이 → 수술 전 대변 마이크로바이옴 → CRC PFS/OS/재발을 formal mediation estimand로 본 코호트는 어디에도 없다. 아시아 CRC 환자 코호트에서 에너지 보정된 식이와 생존을 연결한 논문도 없다. 우리는 FFQ(115명)·수술 전 16S(333명)·enterotype·PFS/OS를 이미 갖고 있어 이 공백을 가장 빨리 메울 수 있는 위치에 있다.

---

## 2. 에너지(칼로리) 보정: 어떻게 할 것인가

### 2-1. 왜 필요한가 (Willett & Stampfer 1986; Willett, Howe, Kushi 1997; Hu 1999; Tomova 2022)
- 총 에너지는 체격·활동량·대사효율과 함께 움직이고, 거의 모든 영양소·식품 섭취량과 강하게 상관한다 → 보정하지 않으면 "많이 먹는 사람"의 효과를 영양소 효과로 오인한다.
- FFQ의 과소/과대 보고 오차는 모든 항목이 공유한다 → 에너지 보정이 이 공통 오차를 상당 부분 제거한다.
- 보정 모델에 따라 추정하는 "질문"이 다르다(Tomova 2022): residual/standard model = **등열량 치환 효과**(다른 열량원 대신 그 영양소를 먹을 때), partition model = **추가 섭취 효과**, density model은 해석이 모호. Letter에서는 어떤 추정량인지 한 문장으로 명시한다.

### 2-2. 우리 데이터의 현재 상태 (반드시 먼저 확인)
| 파일 | 상태 | 조치 |
|---|---|---|
| `1-Input/diet/Dietary_data_2022_115_residual_from_Image06.csv` (115명, 2022년 사용) | 영양소·식품군은 **이미 residual 보정**(209개 중 157개 변수에 음수 존재); `var1`(총열량)은 raw kcal | **재보정 금지.** 다만 원 보정에 로그변환·성별이 들어갔는지 신애선 교수팀(Tung Hoang)에게 확인하고, 가능하면 **raw g/day 원본을 요청**해 아래 2-3을 우리 손으로 다시 수행하는 것이 최선 |
| 540명 임상+FFQ 엑셀 (lab server, git 제외) | raw g/day, kcal | 아래 2-3 절차 적용 |

### 2-3. 권장 절차 (영양소·식품군·식품항목을 구분)
1. **에너지 타당성 필터**: 총열량 <500 또는 >4,000 kcal/일(여), <800 또는 >4,200(남) 제외, 또는 log(에너지) ±3 SD. 제외 수를 Table 1 각주에 보고.
2. **영양소(var2–var94)**: Willett residual method를 **로그 스케일**에서 수행.
   `log(x + ε) ~ log(E)` 회귀 잔차 + `mean(log(x + ε))` → 지수화하면 "평균 에너지 섭취자 기준 섭취량" 단위가 되어 해석이 쉽다. ε은 0이 있는 변수에만 최소 양수값의 절반.
   성별 차이는 (a) Cox 모델에 Sex 공변량, 또는 (b) 성별로 잔차를 따로 구함(Van Blarigan 2020처럼 성별 분위). 115명 규모에서는 (a) 권장.
3. **식품군(group1–17)·식품항목(Food1–106)**: 0 섭취자가 많아 선형 잔차가 잘 맞지 않는다. 주 분석과 민감도 분석을 나눈다.
   - 주: **섭취량 삼분위(또는 비섭취/중간/상위) + 총열량(log kcal)을 Cox 공변량**으로 넣는 standard multivariate model (Van Blarigan 2018, van Zutphen 2021 방식).
   - 민감도 1: `log1p(x)` 잔차(위 2번과 같은 방식).
   - 민감도 2: 밀도(g per 1,000 kcal) + 총열량 공변량(multivariate nutrient-density model).
   - Food46(고사리·고구마줄기·토란대)처럼 저빈도 항목은 주 분석을 이분(비섭취 vs 섭취)으로 두고 위 두 민감도로 방향 일치만 확인한다.
4. **식이 패턴**: 보정된 17개 식품군을 표준화 → PCA + varimax(고유값 >1.5, 적재량 ≥0.2로 이름 붙이기; Meyerhardt 2007, Hoang 2022 방식) → 패턴 점수 삼분위 → Cox에 총열량 공변량 추가. 2–3개 패턴이 나오면 "식물성/전통(채소·해조·콩)" vs "서구/육류" 축이 기대된다.
   - 추가 아이디어(novelty): **RRR(reduced rank regression)** – 반응변수를 Segatella/Leyella/Prevotella CLR로 두고 "마이크로바이옴 지향 식이패턴"을 뽑아 예후와 연결. R `pls` 또는 `rrr` 패키지.
5. **선험적 지수(보완)**: KHEI(한국인 건강식생활지수), aMED, WCRF/AICR 점수, DII(반드시 1,000 kcal당 영양소로 계산한 energy-adjusted DII). 각각 총열량 공변량 포함.
6. **보고 규칙**: Methods에 (i) 보정 모델명, (ii) 로그변환 여부, (iii) 추정량(등열량 치환 효과), (iv) 에너지 제외 기준, (v) 민감도 분석 결과 일치 여부를 적는다. 메타분석 호환을 위해 1 SD당 HR과 삼분위 HR을 함께 보고한다.

### 2-4. 바로 쓸 수 있는 R 스니펫
```r
## residual method on the log scale, mean added back
resid_adj <- function(x, logE, eps = NULL) {
  if (is.null(eps)) eps <- ifelse(any(x == 0, na.rm = TRUE), min(x[x > 0], na.rm = TRUE) / 2, 0)
  y <- log(x + eps); r <- residuals(lm(y ~ logE, na.action = na.exclude))
  exp(r + mean(y, na.rm = TRUE)) - eps          # back to original units
}
dat <- dat %>% filter(var1 >= 500, var1 <= 4000) %>% mutate(logE = log(var1))
dat <- dat %>% mutate(across(c(starts_with("group"), var2:var94), ~ resid_adj(.x, logE), .names = "{.col}_r"))
## standard multivariate model for a food item: tertile + log energy
coxph(Surv(PFS_time, PFS) ~ ntile(Food46, 3) + logE + age + sex + stage34, data = dat)
```
(`5-Reanalysis-2026/2-script/01_diet_prognosis_reanalysis.R`의 Step 1이 같은 논리로 짜여 있다.)

---

## 3. 마이크로바이옴을 매개자로: 분석 설계

### 3-1. 인과 다이어그램(DAG)
- 노출 X: 식이(주: 식이패턴 점수 또는 Food46/채소·해조류 식품군; 보조: 섬유 var7)
- 매개 M: 수술 전 대변 마이크로바이옴 – (i) Segatella/Leyella/Prevotella CLR, (ii) enterotype(이분), (iii) M5 hazard score
- 결과 Y: PFS(주), OS(보조)
- 교란 C: 나이, 성별, 병기(0–2 vs 3–4), 총열량, BMI, 흡연·음주, 신보조치료(식이와 마이크로바이옴 모두에 선행)
- 주의: M5 자체가 X의 하류이자 Y의 예측자 → "M5의 예후 효과가 식이 보정 후 유지되는가"는 독립성 검정이지 매개 검정이 아니다. 둘 다 보고한다.

### 3-2. 분석 단계와 표본
| 단계 | 질문 | 표본 | 방법 |
|---|---|---|---|
| A | X → M | 114–115명 | CLR 값에 선형회귀(공변량 C), Spearman; 다중검정 BH. 조성자료이므로 CLR 또는 ALDEx2 |
| B | M → Y, X → Y | 333명(M→Y), 114명(X→Y) | Cox(병기 보정); X→Y는 2-3의 모델 |
| C | 매개 | 114명, PFS 22건 | 아래 3-3 |
| D | 독립성 | 114명 | Segatella/M5 HR을 식이 보정 전후 비교(이미 handover 결과: 유지됨) |

### 3-3. 매개분석 방법 (권장 순서)
1. **저차원 counterfactual mediation (주 분석)**: R `CMAverse::cmest()` – 매개변수 모델은 선형(CLR), 결과 모델은 Cox 또는 AFT(Weibull), `model = "rb"`(regression-based) + 부트스트랩 CI(≥1,000회). 산출: NDE, NIE, total effect, proportion mediated. Shimomura 2023과 같은 추정량이라 비교·인용이 쉽다. 노출은 연속(1 SD) 또는 이분(상위 vs 하위) 두 가지로.
2. **고차원 탐색 (보조)**: `HIMA::hima_survival()`(또는 `hima(..., mediator.type = "compositional")`)로 70속 CLR을 한꺼번에 스크리닝. 결과는 "탐색적" 꼬리표.
3. **단순 교차검증**: product-of-coefficients(β_XM × β_MY, log HR)와 부트스트랩 – 결과가 1과 같은 방향인지 확인.
4. 해석 원칙: 22건의 이벤트로는 매개 비율의 신뢰구간이 넓다. letter에서는 **점추정 + CI + 방향 일치**를 보여주고, 검정력 한계를 첫 문장에서 인정한다.

### 3-4. 예상 결과 시나리오와 메시지
- 매개 ○: "채소·해조류 식이 → Segatella/Leyella 증가 → 예후 개선"으로 **수술 후 식이 중재의 근거**.
- 매개 ✕(독립): 마이크로바이옴 바이오마커(M5)가 식이와 독립적인 예후인자임을 강화 → 2022 논문의 가치 상승. 어느 쪽이든 letter가 된다(handover 슬라이드 5의 논리 유지).

### 3-5. Letter 형식 제안
- 분량 1,200–1,500단어, Figure 1개(4패널: A 식이→taxa 산점/히트맵, B 식이 삼분위 KM, C 매개분석 forest(NDE/NIE), D M5 HR 식이 보정 전후), Table 1개(모델별 HR), 참고문헌 15–20.
- 후보 저널: *Gut Microbes*(Brief Report), *Microbiome*(Short Report) 또는 *BMC Cancer*(Research note); 외과 쪽 공동저자 선호에 따라 *Annals of Coloproctology*(Letter/Short communication)도 가능.

---

## 4. 강다영 학생 작업 계획 (10주)
| 주차 | 할 일 | 산출물 / 체크포인트 |
|---|---|---|
| 1–2 | 저장소 clone, `run_all.R` 실행으로 2022 논문 재현 확인; Huh 2022·Hoang 2022 정독; 데이터 구조(NO. 키, 115 vs 540, 중복 불가) 숙지 | README 체크리스트 수치 일치 보고 |
| 3 | `02_retaxonomy_prevotella.R`: SILVA 138.2 + GTDB r220으로 ASV 재분류, 잠정 crosswalk와 비교, enterotype 재도출 | crosswalk 표, 일치율 |
| 4–5 | 에너지 보정 파이프라인(2-3), Table 1(식이·임상 기술통계), 식이패턴 PCA | `cox_food46_models.csv`, 패턴 적재량 표 |
| 6 | X→M, M→Y, X→Y Cox 모델(2-3, 3-2 A/B/D) | 모델별 HR 표 |
| 7 | 매개분석(CMAverse + 부트스트랩, HIMA 탐색) | NDE/NIE/PM 표, forest |
| 8 | 그림·표 완성, Methods 초안 | Figure 1, Table 1 |
| 9 | Letter 초안(Introduction 1문단, Results 3문단, Discussion 2문단) | 초안 v1 |
| 10 | 교수 검토 반영, 공동저자(박지원·김민정·신애선 교수팀) 회람 | 투고본 |

**교수가 먼저 결정/확인해야 할 사항**
1. 540명 엑셀의 식이–예후 분석 사용 범위(IRB 범위, 공동연구 합의).
2. 신애선 교수팀에 (a) residual 보정 세부(로그변환·성별), (b) raw g/day 원본 제공 가능 여부 문의.
3. 주 노출을 Food46(단일 식품)으로 할지 식이패턴/식품군(일반화 가능)으로 할지 – 권장: 식이패턴을 주, Food46을 보조.
4. 매개분석 추정량(NDE/NIE on hazard-ratio scale vs survival-time scale) 선택.

## 5. 리스크
- 2022년 FFQ 파일의 보정 세부를 확인하지 못하면 Methods를 정확히 쓸 수 없다 → 2번 확인이 선행 조건.
- 115 vs 114 불일치(1명 대변 미매칭) 기록.
- Prevotella 분할 결과가 잠정값과 다르면 Segatella/Leyella 관련 수치 전부 갱신.
- 이벤트 22건: 매개 효과는 탐색적; 결과 과장 금지.
