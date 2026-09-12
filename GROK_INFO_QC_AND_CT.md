# Квадратичный Чабо (rank 2 / квадратичное поле) и CasselsTatePairing

Grok, 11.09.2026. Только источники и сигнатуры. Не доказательство.

---

## 0. Важная развилка

Готовый QC-код считает **целые точки** \(E(\mathcal O_K)\) или \(E(\mathbb Z)\), не срез
\(\{P\in E(K): x(P)\in\mathbb Q\}\).

Последнее — другая задача: это не \(E(\mathcal O_K)\) и не \(A(\mathbb Q)\) для \(A=\mathrm{Res}_{K/\mathbb Q}E\).
Если коэффициенты Вейерштрасса лежат в \(\mathbb Q\), то \(x\in\mathbb Q\), \(y\in K\) сводится к двум кривым над \(\mathbb Q\) (\(E\) и твист на \(\mathrm{disc}(K)\)). Если коэффициенты в \(K\setminus\mathbb Q\) (как E₁ семейства сечений) — готового солвера нет.

---

## 1. Готовый код QC

### 1.1. То, что ближе всего к «E(K), rank 2, K квадратичное»

**Francesca Bianchi**, Sage, мнимое квадратичное, E задана над \(\mathbb Q\),
\(\mathrm{rank}\,E(\mathbb Q)=1\), \(\mathrm{rank}\,E(L)=2\):

- https://github.com/bianchifrancesca/QC_elliptic_imaginary_quadratic_rank_2
- Главная функция: `quad_chab_ell_im_quad(E, p, n, double_root_prec, L, ...)`
  в `quad_chab_ell_im_quad.sage`
- Статья: Balakrishnan–Besser–Bianchi–Müller, *Explicit quadratic Chabauty over number fields*,
  Israel J. Math. 243 (2021) = arXiv:1910.04653. §6.1–6.2: эллиптические кривые над вещественным и мнимым квадратичным.
- Ограничения кода: плохие простые с Tamagawa >1 должны быть **инертны или разветвлены** в L;
  p — хорошая редукция и **расщеплена** в L.
- Пример из docstring: `E = EllipticCurve("91a1")`, `K.<a> = QuadraticField(-1)`, `p=5`.
- Выход: целые точки \(E(\mathbb Z)\) и \(E(\mathcal O_L)\), плюс лишние p-адические.

**Aashraya Jha** (Sage 10.2, порт Bianchi; E/K мнимое квадратичное, rank 2):

- https://github.com/AashrayaJha/QC_ECIm
- Статья: arXiv:2311.01691
- Считает p-адические **целые** точки; решето, если кручение тривиально.
- Главные файлы: `Quadratic_Chabauty_for_ECs_2/quad_chab_ell_im_quad.sage`, `Heights_for_im_quad.sage`.

**Bianchi над \(\mathbb Q\), rank ≤1** (не rank 2 над K):

- https://github.com/bianchifrancesca/quadratic_chabauty
- `quadratic_chabauty_elliptic.sage`: `quadratic_chabauty_rank_0`, `quadratic_chabauty_rank_1`
- Статья: ANT 14 (2020) = arXiv:1904.04622. Считает \(X(\mathbb Z)\), не \(E(K)\).

### 1.2. Род 2 / биэллиптические (не E/K)

- Bianchi–Padurariu: https://github.com/bianchifrancesca/QC_bielliptic
  `qc_g2_bielliptic.sage` — род 2 над \(\mathbb Q\), ранг 2, уравнение \(y^2=a_6x^6+\cdots+a_0\), два эллиптических фактора ранга 1.
- Gajović–Müller LinQC: https://github.com/steffenmueller/LinQC
  `LinQCRQF.sage` — род 2 над **вещественным** квадратичным, целые точки. Нужны Sage+Magma и
  https://github.com/jbalakrishnan/AWS
- QCMod (Magma, модулярные кривые над \(\mathbb Q\), rank=genus): https://github.com/steffenmueller/QCMod
- QC над \(\mathbb Q(\zeta_3)\), род 3, rank 6: arXiv:2501.07833 (Balakrishnan–Betts–Hast–…–Müller). Код к статье, не универсальный солвер E/K.

### 1.3. Magma встроенный `Chabauty`

Это **эллиптический** Чабо (ранг 1), не квадратичный:

```
Chabauty(mwmap, EtoP1)
```

Handbook: elliptic curve Chabauty, Bru04. Для rank E(K)=2 **не** применим.

---

## 2. CasselsTatePairing в Magma: сигнатуры

Handbook (Curves over Q / number fields), V2.15+:

```
CasselsTatePairing(C, D) : CrvHyp, CrvHyp -> RngIntElt
```

- C, D — гиперэллиптические **y² = квартика** (степень 4), 2-накрытия **одной** эллиптической кривой.
- База: \(\mathbb Q\), **числовое поле**, или \(\mathbb F(t)\) (char ≠ 2).
- Обе кривые должны быть локально разрешимы во всех пополнениях, иначе спаривание не определено.
- Значение в \(\mathbb Z/2\mathbb Z\), возвращается как 0 или 1.
- Типичный источник C, D: `TwoDescent(E)`.
- Verbose: `SetVerbose("CasselsTate", 1)` или `2`.

Вторая сигнатура — **только над \(\mathbb Q\)**:

```
CasselsTatePairing(C, D) : Crv, CrvHyp -> RngIntElt
```

C — 4-накрытие (пересечение двух квадрик в \(\mathbb P^3\), из `FourDescent`), D — 2-накрытие y²=квартика.

`TwoPowerIsogenyDescentRankBound` — только `CrvEll[FldRat]`, не квадратичное поле.

Связанное (не то же самое):

```
TwoDescent(E) : CrvEll -> [CrvHyp], [Map], Map     // над Q и над числовым полем
TwoSelmerGroup(E)
TwoCover(e)                                        // одно 2-накрытие из элемента алгебры
```

Над числовым полем `TwoDescent` тянет группу классов/единицы. Ускорение:

```
SetClassGroupBounds("GRH");
```

---

## 3. Примеры Magma

### 3.1. Гарантированно < 1 с (над \(\mathbb Q\))

Handbook example, 0.27 + 0.13 с:

```
E := EllipticCurve("571a1");
covers := TwoDescent(E);
CasselsTatePairing(covers[1], covers[2]);  // -> 1
```

Или сразу Sha:

```
MordellWeilShaInformation(EllipticCurve("571a1") : ShaInfo);
```

(в handbook ~0.84 с; печатает матрицу CT на Sel₂.)

### 3.2. Над квадратичным полем, цель ≤ 60 с

Спаривание само по себе дешёвое; лимит съедает `TwoDescent`. GRH обязателен.
Минимальный payload (класс числа 1, 2-кручение, маленькие коэффициенты):

```
SetClassGroupBounds("GRH");
K<r> := QuadraticField(-3);
E := EllipticCurve([K| 0, 0, 0, 1, 0]);  // y^2 = x^3 + x
covers := TwoDescent(E);
print "ncovers =", #covers;
if #covers ge 2 then
  time print CasselsTatePairing(covers[1], covers[2]);
end if;
```

Если `TwoDescent` не уложится: не строить все модели, а

```
SetClassGroupBounds("GRH");
K<r> := QuadraticField(-3);
E := EllipticCurve([K| 0, 0, 0, 1, 0]);
G, m := TwoSelmerGroup(E);
print G;
```

`TwoSelmerGroup` обычно быстрее полного списка приведённых квартик; CT требует именно квартики (`covers[i]`).

Кривая из MathOverflow (может быть тяжелее 60 с):

```
SetClassGroupBounds("GRH");
K<b> := QuadraticField(-37);
A := EllipticCurve([K| 0, 0, 0, 17, 0]);
covers := TwoDescent(A);
```

Для семейства сечений (E₁/ℚ(√d) с большими a₄, a₆) `TwoDescent` на онлайн-калькуляторе с лимитом 60 с, скорее всего, не закончится — как уже было с `MordellWeilShaInformation(E1)`. Практичнее CT на **якобиане рода 2 над ℚ** (Fisher–Yan, github.com/TomAFisher/Cassels-Tate-pairing--Genus-2: `CasselsTatePairingGenus2` после attach; это не встроенный `CasselsTatePairing`).

---

## 4. Что из этого брать для сечений с w = +1

Готовый `quad_chab_ell_im_quad` **не** подходит к E₁ семейства:

- E₁ задана над K, не спуск с \(\mathbb Q\) (`descend_to(QQ)` пуст);
- поля вещественные (√3, √10, √42, …), код Bianchi/Jha — для **мнимых**;
- плохие простые E₁ не обязаны быть инертны;
- цель — x∈ℚ, не \(\mathcal O_K\)-целые точки.

Ближний обход, если когда-нибудь E опустится на \(\mathbb Q\): rank E(K)=rank E(ℚ)+rank E^{(d)}(ℚ); при ранге 2 на одном из двух факторов над ℚ — `quadratic_chabauty_rank_1` не хватит, нужен QC на той кривой над \(\mathbb Q\) (rank 2) — готового эллиптического QC rank 2 над \(\mathbb Q\) в репозиториях Bianchi нет (там rank ≤1 над \(\mathbb Q\)).

Конец.
