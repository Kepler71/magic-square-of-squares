# Проверка новизны по литературе, 2026-09-12

Задание координатора: по шести пунктам вынести вердикт **известно / частично известно / упоминаний не нашёл**,
со ссылками и короткими цитатами. Цель — честная подготовка письма к A. Auel (Dartmouth, arXiv:2609.09351).

Статусы цитат: **цитата** — дословно из источника; **вывод** — моя выкладка; **не найдено** — честно не найдено
(отсутствие результата поиска, а не доказательство отсутствия).

Что реально прочитано целиком в этой проверке (не по аннотациям):
- A. Bremner, *On squares of squares*, Acta Arith. 88 (1999), 289–297 — полный текст (PDF с multimagie.com/Bremner1.pdf
  и pdfs.semanticscholar.org), 9 страниц, прочитан постранично.
- A. Bremner, *On squares of squares II*, Acta Arith. 99 (2001), 289–308 — полный текст (multimagie.com/Bremner2.pdf).
- C. Boyer, *A search for 3x3 magic squares having more than six square integers…* (2004), multimagie.com/Search.pdf — целиком.
- L. Morgenstern, страницы на multimagie.com: Morgenstern08 (окт. 2010), Morgenstern11 (нояб. 2011), Morgenstern13, Morgenstern25,
  и PDF *3x3 Magic Square of Squares Properties* (июль 2015), multimagie.com/MorgensternMssProperties.pdf.
- C. Boyer, страница *Latest research on the «3x3 magic square of squares» problem* — целиком (2006–2020).
- P. Pierrat, F. Thiriet, P. Zimmermann, *Magic squares of squares* (2015), members.loria.fr/PZimmermann/papers/squares.pdf — целиком.
- A. Auel, B. Singer, arXiv:2609.09351 — исходник TeX с arXiv (e-print), грепы по родам, Чабо, спуску.
- E. González-Jiménez, arXiv:1311.5759, §1–2; M. Stoll, arXiv:1711.00500, введение и §1.
- Локальная проверка Sage 10.x / PARI 2.17.4 на наличие спаривания Касселса–Тейта.

---

## 1. Лемма о двух полных парах

**Наше утверждение.** Всякий магический 3×3 с ≥ 7 квадратными клетками содержит не менее двух полных пар
противоположных клеток; следствие — семейства с двумя полными парами (G-семейства) покрывают задачу о семи квадратах
полностью по постановке.

### Вердикт: **ИЗВЕСТНО** (и сама лемма, и следствие).

**(а) Сама лемма — у Бойера, 2004, дословно.**
C. Boyer, *A search for 3x3 magic squares having more than six square integers among their nine distinct integers*,
Draft v2, 16.09.2004, http://www.multimagie.com/Search.pdf, раздел «Some words about the method used», с. 3:

> «A line going through the central cell C, and having two square integers around the central cell, is an integer
> solution of the equation: x² + y² = 2C. … **All the 7.x and 8.x configurations need two, three or four such lines
> through the centre, meaning at least (as a strict minimum) two, three or four solutions of the above equation.**»

Это ровно наша лемма, сформулированная для всех восьми конфигураций 7.I–7.VIII (там же, fig. 5 — полный список
восьми конфигураций семи клеток с точностью до симметрий D₈) и трёх конфигураций 8.I–8.III.
Бойер идёт дальше и выписывает счётные формулы (D1/D2): при квадратном центре c² с n различными простыми 4k+1
число решений x²+y²=2c², x<y, равно (3ⁿ−1)/2; при неквадратном центре C — 2^(n−1).

**(б) Структурная основа — у Моргенштерна, 2015, дословно.**
L. Morgenstern, *3x3 Magic Square of Squares Properties*, July 2015, http://www.multimagie.com/MorgensternMssProperties.pdf,
Lemma 9:

> «**Eight APs exist in a 3x3 magic square of squares. Four of the APs run through the center and cover all nine entries.**
> x−y, x, x+y (step y); x−z, x, x+z (step z); x−y−z, x, x+y+z (step y+z); x−y+z, x, x+y−z (step y−z).
> Four more APs are on the pandiagonals.»

Четыре AP через центр — это наши четыре полные пары; они попарно не пересекаются и покрывают все девять клеток.
Отсюда «две неквадратные клетки портят не более двух пар» — одна строка. (Работа явно указывает, что она переписывает
элементарно результаты L. W. Rabern, *Properties of magic squares of squares*, Rose-Hulman Undergrad. Math. J. 4 (2003),
https://scholar.rose-hulman.edu/rhumj/vol4/iss1/3/.)

**(в) Следствие («две полные пары покрывают задачу о семи квадратах») — у Моргенштерна, 2010, дословно и с полнотой.**
L. Morgenstern, *Search method to find 3x3 magic squares with 7 distinct square entries*, October 2010,
http://www.multimagie.com/English/Morgenstern08.htm:

> «This method is **complete, non-redundant** … **Any 7-square solution can be formed from one of the following six
> configurations of 5-square solutions.** Configurations 1,2,3 correspond to 7-square solutions where the central entry
> is a square. Configurations 4,5,6 correspond to 7-square solutions where the central entry is not necessarily a square.»
> «Each of the configurations of 5 squares consists of **two 3-square arithmetic progressions having one entry in common**.
> 1: A²+I² = C²+G² = 2E²;  2: H²+B² = F²+D² = 2E²;  3: A²+I² = H²+B² = 2E²;
> 4: 2G²−F² = 2I²−D² = B²;  5: 2A²−F² = 2C²−D² = H²;  6: 2G²−B² = 2A²−H² = F².»

Конфигурации 1, 2, 3 — это буквально наши G2 (две угловые пары + центр), G1 (две рёберные пары + центр),
G3 (угловая + рёберная + центр). Конфигурации 4–6 — случай неквадратного центра.
Параметризация у него та же, что наша Сегре-параметризация: X=2mn, Y=m²−n², Z=m²+n², AP = ((Y−X)², Z², (Y+X)²),
и «two progressions having a common entry by scaling each progression by the common entry of the other».

**Что здесь всё же наше.** (i) Формулировка как отдельной леммы с одностроковым доказательством — да, но это не результат.
(ii) Наш вклад не в лемме, а в том, что́ с G-семейством дальше делается (кривая рода 5 и арифметика на ней).
Заявлять лемму как новую **нельзя**.

---

## 2. Структура семейства: b = 1+2s−s², h = 1−2s−s², n = 1+s²

**Наше утверждение.** C²−A² = −8s(s−1)(s+1)(s²+1)²; поле сечения всегда k_s = ℚ(√(−8s(s²−1)));
база семейства — кривая рода 1 (твист y²=x³−x); j(E₁) ∈ ℚ(s); полное 2-кручение над k_s всегда.

### Вердикт: **ЧАСТИЧНО ИЗВЕСТНО**. Параметризация и кубика s³−s — классика и есть в печати, причём с теми же
### буквами; специфические арифметические следствия (j ∈ ℚ(s), полное 2-кручение, формула для C²−A²) я не нашёл.

**(а) Параметризация (b,h,n) — классическая, три независимых источника.**

- L. Morgenstern, *3x3 MSS Properties* (2015), Lemma 7 (цитата):
  > «A primitive AP has the formula a = 2mn − m² + n², b = 2mn + m² − n², c = m² + n², with m and n coprime, one odd, one even.»
  При s = n/m это в точности (−h, b, n) с точностью до знака и множителя m².
- P. Pierrat, F. Thiriet, P. Zimmermann, *Magic squares of squares* (2015), Theorem 1 (цитата):
  > «Let p be a square-free divisor of A, p ≡ 1 mod 4. Write A = pA′ … Then write b = 4mn(m²−n²),
  > x = √(A²−p²b), y = √(A²+p²b).» и в доказательстве: «A²−p²b and A²+p²b are perfect squares, respectively of
  > **x = p(m²−2mn−n²)** and of **y = p(m²+2mn−n²)**.»
- Auel–Singer, arXiv:2609.09351, §3.1 (цитата): «with the relation B² + C² = 2A², i.e., **the conic parametrizing squares
  in arithmetic progression**. This is a smooth conic with a rational point [1:1:1], so it is a rational curve on V.»

**(б) Нормированный шаг s(s²−1)/(s²+1)² — дословно у Моргенштерна, 2011.**
L. Morgenstern, *Magic Hourglass Anti-Closure*, November 2011, http://www.multimagie.com/English/Morgenstern11.htm:

> «**The Magic Hourglass Set is the set of positive rational numbers generated by the expression p(1−p²)/(1+p²)²
> for rational p, where 0 < p < 1.** Proving that this set is anti-closed under addition is equivalent to proving that
> there is no magic hourglass. Finding a counterexample to anti-closure is equivalent to finding a magic hourglass.
> Find rational numbers p, q, r such that p(1−p²)/(1+p²)² + q(1−q²)/(1+q²)² = r(1−r²)/(1+r²)².»

**Вывод (проверено выкладкой):** шаг нашей пары равен (b²−h²)/2 = −4s(s²−1), центр (1+s²)², то есть
нормированный шаг = −4·s(s²−1)/(s²+1)² — ровно образующая «магического множества песочных часов» Моргенштерна.
Значит наш параметр s и его роль — это его p, а условие «песочные часы» — его уравнение анти-замкнутости.

**(в) Кубика s³−s и кривая y²=x³−x — есть и у Бремнера (1999), и у Ауэля–Сингера (2026), и восходит к Робертсону (1996).**

- Бойер (M.I. article, начало, цитата): «**John P. Robertson [51] showed that the problem is equivalent to other
  mathematical problems on arithmetic progressions, on Pythagorean right triangles, on congruent numbers and
  elliptic curves y² = x³ − n²x.**» (J. P. Robertson, *Magic squares of squares*, Math. Mag. 69 (1996), 289–293.)
- **Bremner 1999, последняя формула статьи (с. 296) — это буквально наши b, h, n и наша кубика.** Цитата:
  > «or indeed, the family
  > [ (µ²+1)²  −(µ²+2µ−1)²  4(µ³−µ) / −(µ²−2µ−1)²  0  (µ²−2µ−1)² / −4(µ³−µ)  (µ²+2µ−1)²  −(µ²+1)² ]
  > **over the field ℚ(i, √(µ³−µ)).**»
  Здесь µ²+1 = n, µ²+2µ−1 = −h, µ²−2µ−1 = −b, 4(µ³−µ) = −(шаг). То есть Бремнер уже в 1999 г. писал
  наше семейство и указал его поле как ℚ(i, √(µ³−µ)) — а µ³−µ = s(s²−1) — свободная часть нашего k_s
  (наше k_s = ℚ(√(−8s(s²−1))) отличается на √2 и знак; это тот же congruent-number кубический множитель).
  Отличие: у Бремнера центр равен нулю (M = 0), у нас центр ненулевой.
- Auel–Singer, arXiv:2609.09351, Prop. 3.3, 3.4 (цитаты):
  > «The M = 0 hyperplane section decomposes, over ℚ(i), into the union of **16 smooth elliptic curves all isomorphic
  > to y² = x³ − x**.»
  > «For any squarefree integer d, the elliptic curves y² = x³ − d²x are in the same twist family … has positive rank
  > over ℚ if and only if d is a congruent number. Hence, each irreducible component of the M = 0 locus over ℚ(i) has
  > infinitely many rational points over ℚ(i, √d) whenever d is a congruent number.»

**(г) Чего не нашёл (не найдено).**
- Формулы C² − A² = −8·s(s−1)(s+1)·(s²+1)² в печати.
- Утверждения «поле сечения всегда ℚ(√(−8s(s²−1)))», то есть что база *нашего* однопараметрического семейства
  сечений — кривая рода 1, а не прямая.
- Символьного j-инварианта j(E₁) ∈ ℚ(s) и вывода «E₁ — всюду твист кривой над ℚ(s)».
- Утверждения о полном 2-кручении над k_s во всём семействе.
- Статистики знака функционального уравнения по семейству (≈ 50/50).

**Честный итог по п. 2:** «база — твист y²=x³−x, поле связано с s³−s» — это переоткрытие того, что уже стоит
у Бремнера (1999) и у Ауэля–Сингера (2026). Новое здесь — только уровень «арифметики слоя»: j ∈ ℚ(s), полное
2-кручение, следствие о неприменимости квадратичного Шаботи ко всему семейству. Это надо подавать именно так.

---

## 3. Роды: g = 1 + 2^{k−1}(k−2), то есть 0, 1, 5, 17, 49

### Вердикт: **ЧАСТИЧНО ИЗВЕСТНО** — все пять чисел уже напечатаны (Ауэль–Сингер, 08.09.2026),
### и у Бремнера есть один явный род; замкнутой формулы по k в этом контексте я не нашёл.

**(а) Бремнер 2001: «высокий род» — но не только.** Цитата (с. 290):

> «This is done by parametrizations in one variable, so that asking for a seventh entry to be square in these examples
> involves finding rational points on **hyperelliptic curves (in general of high genus)** of type f(t) = □.»

Однако в той же статье он выписывает **явный род** (с. 295, Remark к конфигурации II):

> «For m₁₁ or m₃₃ to be square, the condition is that **the curve of genus 3** given by
> 1 + 4x + 8x² − 28x³ − 2x⁴ − 84x⁵ + 72x⁶ + 108x⁷ + 81x⁸ = □ have rational points.»

Плюс он приводит степени параметризаций всех 16 конфигураций (с. 296): «respective degrees 12, 8, 12, 20, 12, 12, 12,
20, 20, 12, 20, 20, 8, 12, 20, 12». Для гиперэллиптической кривой f(t)=□ с f степени d род = ⌈d/2⌉ − 1, то есть
из его таблицы читаются роды 5 (d=12), 3 (d=8), 9 (d=20). **Формулировка «у Бремнера сказано лишь: кривые вообще
высокого рода» — завышение нашей новизны.** Точнее: он даёт «в общем случае высокого рода» как общее замечание,
но в разобранном примере называет род 3 явно, и его степени параметризаций дают роды напрямую.

**(б) Ауэль–Сингер 2026: все пять чисел 0, 1, 5, 17, 49 в печати.**
arXiv:2609.09351 (исходник TeX), Prop. 3.1 (таблица): компоненты вырожденного локуса Z₃/Z₅/Z₇ имеют
**Genus 0 / 5 / 49** при степенях 2 / 8 / 32. Prop. 3.6 (таблица): **Genus 1 / 5 / 17** при степенях 4 / 8 / 16.
Prop. 3.7: «**degree 8, genus 5 curves**», «**degree 16, genus 17 curves**», «**degree 32, genus 49 curves**».
§4.1: «Their irreducible components are **smooth curves of genus 17 and degree 16, or genus 5 and degree 8**.»
Remark 3.2 (цитата): «The fact that the components of Z₅ and Z₇ are of high genus is consistent, via Faltings theorem,
with older results to the effect that there are no integer 3×3 magic squares of squares with 5 or 7 distinct entries.»

Совпадение не случайное: и у них, и у нас это (ℤ/2)^k-накрытия рациональной базы степени 2^k.

**(в) Формула как таковая.** Не найдено в контексте магических квадратов. Но она — одна строка Римана–Гурвица:
для (ℤ/2)^k-накрытия P¹, заданного k независимыми квадратичными условиями u_i² = F_i (deg F_i = 2, 2k точек ветвления,
инерция ℤ/2), 2g − 2 = 2^k·(−2) + 2k·2^{k−1} ⇒ **g = 1 + 2^{k−1}(k−2)**. Последовательность 0, 1, 5, 17, 49, 129, …
это OEIS **A000337**, a(n) = (n−1)2ⁿ + 1 (там же комментарий: «Genus of graph of n-cube = a(n−3)»).
То есть формула — стандартное упражнение, а не результат; ценность только в том, что она привязана к «k дополнительных
квадратных клеток».

**Честный итог по п. 3:** род 5 для седьмой клетки — **не новость** после 08.09.2026; в письме это надо признать явно
и спросить, совпадают ли наши кривые рода 5 с их компонентами.

---

## 4. Расщепление якобиана рода 5 (1+1+1+2) и использование фактора ранга 0

**Наше утверждение.** Для C: u₁²=F₁(t), u₂²=F₂(t), u₃²=F₃(t) (F_i квадратичные) J_C ~ J₁₂×J₁₃×J₂₃×J_H,
где C_ij: v²=F_iF_j эллиптические, H: Y²=F₁F₂F₃ рода 2; фактор ранга 0 даёт полное перечисление рациональных точек
без Шаботи.

### Вердикт: **ЧАСТИЧНО ИЗВЕСТНО** — расщепление классическое и уже применялось к родственным кривым рода 5,
### а трюк «фактор ранга 0 закрывает семейство» применён к ЭТОЙ задаче лично Бремнером в 2001 г.

**(а) Трюк «фактор ранга 0» — у Бремнера, 2001, в самой этой задаче. Дословно (с. 295):**

> «From (10), the square at (9) has m₂₂ a perfect square provided the curve (1+X)(1+9X)(1−2X+9X²) = □ has rational
> points with X = µ². **But this elliptic curve of conductor 48 has rank 0, forcing µ = 0 and a consequent trivial square.**»

Это в точности схема «эллиптический фактор ранга 0 ⇒ нет нетривиальных точек ⇒ бесконечное семейство закрыто»,
применённая к магическому квадрату квадратов.

**(б) Расщепление якобиана кривой высшего рода на эллиптические факторы — у него же, там же:**

> «**The curve contains in its Jacobian the elliptic curve U⁴ + 4U³ − 4U² − 64U − 32 = □ (observe the transformation
> U = 3x + 1/x), of conductor 1104 and rank 1**, with (−2, 8) a generator for the rational points. But computing the
> “small” rational points on this curve led to no non-trivial rational points on the curve of genus 3.»

То есть Бремнер уже делал: род 3 → эллиптический фактор в якобиане → попытка перечислить точки. И честно писал, что
не смог: «it seems likely that x = 0, ∞ are the only such, though we are unable to show this».

**(в) Род 5 именно такого типа — отдельная литература.**
- E. González-Jiménez, *Covering techniques and rational points on some genus 5 curves*, arXiv:1311.5759 (Acta Arith.),
  §1 (цитата): «Finally, we have the case C : {aX₀²+bX₁²=X₂², cX₀²+dX₁²=X₃², eX₀²+fX₁²=X₄²}. **This curve is generically
  of genus 5 and there are not known algorithms to compute C(K).** In this paper, our purpose is to give an algorithm
  to compute (under some hypotheses) C(K).» Его алгоритм — covering collections + elliptic curve Chabauty;
  попутно там же появляется **род 17** для (ℤ/2)²-накрытия D.
- M. Stoll, *Diagonal genus 5 curves, elliptic curves over ℚ(t), and rational diophantine quintuples*, arXiv:1711.00500
  (Acta Arith. 190 (2019)), введение (цитата): «Some important features that we use are **the splitting of the Jacobian J
  of C as a product of five elliptic curves over ℚ, up to isogeny**, and the fact that J has a large rational 2-torsion
  subgroup, which allows us to consider many étale double coverings of C. The Prym varieties of these coverings are
  isogenous (over ℚ̄) to a product of four elliptic curves defined (in general) over a biquadratic number field.
  This allows us to set up various ways of applying ‘Elliptic Curve Chabauty’ [Bru03] to our situation.»
  (В «диагональном» случае, когда все F_i чётные, наш род-2 фактор H дополнительно расщепляется — отсюда у Штолля
  пять эллиптических кривых вместо 1+1+1+2.)
- Также E. González-Jiménez, X. Xarles (цитируемые обоими) — семейства кривых рода 5 тем же методом.

**(г) Чего не нашёл (не найдено).** Применения разложения 1+1+1+2 (или пятерного) именно к кривым, отвечающим
седьмой/восьмой клетке магического квадрата, и вообще любого использования elliptic curve Chabauty, спуска или
Селмеровых групп в литературе о магическом квадрате квадратов. У Ауэля–Сингера ни «Chabauty», ни «Selmer»,
ни «descent», ни рангов Морделла–Вейля кривых нет (грепы по исходнику TeX — пусто).

**Честный итог по п. 4:** метод классический, применялся к похожим кривым рода 5, и его зародыш есть у Бремнера
в этой самой задаче. Новое — конкретное применение в нашей конфигурации (какая именно кривая, какие F_i, какой фактор
оказывается ранга 0). Подавать как «новый метод» нельзя; как «новое применение» — можно.

---

## 5. Исключение бесконечных семейств отношений

**Наше утверждение.** Ни одна рёберная пара магического квадрата из девяти различных квадратов не может иметь
отношение (17:7:13), (23:7:17), (71:49:61) — при любом масштабе.

### Вердикт: **УПОМИНАНИЙ НЕ НАШЁЛ** в этой форме. Это, по-видимому, единственный из шести пунктов,
### где заявка на новизну выдерживает проверку.

**(а) Все известные результаты о пустоте — это границы по центру или по шагу, то есть конечные области.**
- D. Buell, *A search for a magic hourglass* (1998/1999), multimagie.com/Buell.pdf. Цитата по Бремнеру (2001, с. 291):
  «Duncan Buell [5] has shown by careful search that there is no seven-square magic square corresponding to the
  “hour-glass” configuration in which **the central element of the square is less than 25·10²⁴**.»
- C. Boyer (2004, Search.pdf): перебор центров по типам разложения на простые 4k+1, все границы конечные
  (до ~2,5·10²⁷ в самом дальнем случае).
- L. Morgenstern (2008, 2013, 2014): «All primitive APs up to **d = 6.4×10²²**», позже «up to **d = 6.0×10²³**,
  the primitive step values being 4mn(m²−n²), with m,n coprime, one odd, one even, and **n < m < 2²⁴**».
  Здесь ограничение именно на генераторы (m,n) — то есть на *отношение* — но сверху и конечное.
- L. Rabern (2007), R. Rathbun (2006, 2010, 2017), L. Pech (2006), Pierrat–Thiriet–Zimmermann (2015) — то же самое:
  конечные диапазоны, либо решения по модулю 2^N (до 2^59).

**(б) Исключения целых параметрических семейств в литературе есть, но другого типа — «семейство параметризовано,
покажем, что в нём нет решений», а не «фиксировано отношение одной пары, покажем, что нет решений при всех масштабах».**
- Бремнер 2001 (см. п. 4а): кривая кондуктора 48 ранга 0 ⇒ µ = 0 ⇒ тривиальный квадрат. Это исключение
  бесконечного однопараметрического семейства через ранг 0. **Ближайший по духу прецедент.**
- Бойер (страница Latest research, 22.01.2010, о квадрате Ф. Рубина): «It is known that unfortunately **this family
  [Lucas] can't produce a 3x3 magic square of squares**.»
- Бойер о семействе Ant King (27.05.2011): «Because its magic sum is a squared integer …, **this parametric solution
  unfortunately can't produce magic squares of squares** (needing a magic sum equaling three times a squared integer).»
- Бойер о семействе A. Wesołowski (20–22.12.2019): «**When n is an even positive number, we can prove that the cells
  in red can't be squared integers.**» — исключение половины семейства.
  (Замечание: это семейство
  [(xy−z)², (xz+y)², x²+y²z² / (x+yz)², (x²+y²)(z²+1)/2, (xz−y)² / x²z²+y², (yz−x)², (xy+z)²]
  совпадает по форме с нашим G1 — те же три «свободные» клетки F₀ = m²+n²t², F₄ = (m²+n²)(t²+1)/2, F₈ = n²+m²t².
  Это стоит проверить отдельно: возможно, наша G1 и его семейство — один и тот же объект.)

**(в) Отдельно.** Утверждений вида «данное фиксированное примитивное AP (b², n², h²) не достраивается до магического
квадрата ни при каком масштабирующем множителе» я не нашёл ни у Бойера, ни у Моргенштерна, ни у Бремнера,
ни у Ауэля–Сингера, ни в arXiv (проверено поиском по формулировкам «fixed arithmetic progression», «all scalings»,
«infinitely many scalings excluded»).

**Честный итог по п. 5:** формулировка «закрыто бесконечное семейство масштабов фиксированного отношения»
в литературе не встречается; это наш тип результата. Но и вес его надо назвать честно: три отношения из бесконечного
множества, и только для двух рёберных пар.

---

## 6. Открытая реализация спаривания Касселса–Тейта над квадратичным полем с полным 2-кручением

### Вердикт: **УПОМИНАНИЙ НЕ НАШЁЛ** — открытых реализаций нет; вывод проекта подтверждается, в том числе локальной проверкой.

**(а) PARI/GP 2.17.4 (проверено локально, `gp -q`).**
- `ellrank(E,{effort},{points})`: «**if E is an elliptic curve over ℚ**, attempts to compute the Mordell-Weil group…»
- `ell2cover(E)`: «**if E is an elliptic curve over ℚ**, returns a basis of the set of everywhere locally soluble 2-covers…»
- `ellrankinit(E)`: «**if E is an elliptic curve over ℚ**, initialize data for further calls to ellrank.»
- Идентификатора `ellcassels` (или любого с `cassels`) нет.
Спаривание Касселса в PARI используется внутри `ellrank` — но только над ℚ; см. документацию Sage
(`sage/schemes/elliptic_curves/ell_rational_field.py`, `rank_bound`): «The following is the curve 571a1, which has rank 0,
but Sha of order 4, yet **pari, using the Cassels pairing** is able to show that the rank is 0.»

**(б) SageMath (проверено грепом по установленным исходникам).** Никакой функции спаривания Касселса–Тейта нет;
единственные вхождения слова «Cassels» — в комментариях (ссылка на PARI выше и ссылка на учебник Касселса в `hom_velusqrt.py`).

**(в) Magma.** Реализация S. Donnelly, `CasselsTatePairing`. Release notes Magma V2.15:
«An efficient implementation of the Cassels-Tate pairing between 2-coverings is provided, using **an algorithm due to
Donnelly** which is much cleaner than the version previously implemented.» Страница Acknowledgements: S. Donnelly
(2005–2017) — «many routines for elliptic curves **over ℚ and number fields**, including descent methods and
Cassels–Tate pairings». То есть над числовым полем это, по-видимому, работает — но Magma закрытая и платная.

**(г) Академические реализации — все в Magma, публичных репозиториев не найдено.**
- H. Shukla, M. Stoll, arXiv:2302.01640 (над произвольным числовым полем) — кода в открытом доступе не нашёл.
- T. Fisher, arXiv:2208.14977 (бинарные квартики, без решения коник) — программы Фишера распространяются как Magma-файлы
  с его страницы dpmms.cam.ac.uk/~taf1000 (страница отдаёт 403 из этой среды; проверить не удалось).
- J. Yan, arXiv:2109.08258; T. Fisher, J. Yan, arXiv:2306.06011 (род 2 над числовым полем) — Magma, репозиторий не найден.
- T. Fisher, R. Newton (3-Selmer), van Beek–Fisher (3-изогении) — Magma.
- Поиск по GitHub/arXiv на «Cassels–Tate pairing implementation / code» не дал ни одного открытого репозитория.

**Честный итог по п. 6:** наша реализация `descent/ctp.sage` (Sage, над ℚ(√D), полное 2-кручение), похоже, —
единственная открытая. Это **инструментальная**, а не математическая новизна, и подавать её надо именно так.

---

## Итоговая таблица

| № | Утверждение | Вердикт | Кто и где |
|---|---|---|---|
| 1 | Лемма о ≥ 2 полных парах при ≥ 7 квадратах | **не новое** | Boyer, Search.pdf (2004), §«method»; Morgenstern, MssProperties (2015), Lemma 9 |
| 1′ | Следствие: две полные пары покрывают задачу о 7 квадратах | **не новое** | Morgenstern, Morgenstern08.htm (окт. 2010), «complete, non-redundant», конфигурации 1–6 |
| 2a | Параметризация b, h, n конического сечения | **не новое** | Morgenstern Lemma 7; Pierrat–Thiriet–Zimmermann Thm 1; Auel–Singer §3.1 |
| 2b | Нормированный шаг s(s²−1)/(s²+1)² | **не новое** | Morgenstern, Magic Hourglass Anti-Closure (нояб. 2011) |
| 2c | База семейства — твист y²=x³−x, поле через s³−s | **не новое** | Bremner 1999, с. 296 (ℚ(i,√(µ³−µ)) с теми же b,h,n); Auel–Singer Prop. 3.3–3.4; Robertson 1996 |
| 2d | j(E₁) ∈ ℚ(s); полное 2-кручение над k_s всегда; C²−A² = −8s(s²−1)(s²+1)² | **под вопросом** (упоминаний не нашёл) | — |
| 3 | Роды 0, 1, 5, 17, 49 | **не новое** (с 08.09.2026) | Auel–Singer, Prop. 3.1, 3.6, 3.7, §4.1; отдельно Bremner 2001 — явный род 3 и степени 8/12/20 |
| 3′ | Замкнутая формула g = 1 + 2^{k−1}(k−2) в этом контексте | **под вопросом** — сама формула тривиальна (Риман–Гурвиц; OEIS A000337) | — |
| 4 | Расщепление якобиана рода 5 на 1+1+1+2 | **не новое** | классика (ℤ/2)³-накрытий; Stoll arXiv:1711.00500 (5 эллиптических факторов), González-Jiménez arXiv:1311.5759 |
| 4′ | Фактор ранга 0 ⇒ полное перечисление точек, в задаче о магическом квадрате | **не новое** | Bremner 2001, с. 295 (кондуктор 48, ранг 0, «forcing µ = 0») |
| 4″ | Применение этого к кривой рода 5 седьмой клетки | **под вопросом** (упоминаний не нашёл) | — |
| 5 | Исключение отношения пары при ВСЕХ масштабах: (17:7:13), (23:7:17), (71:49:61) | **новое** (упоминаний не нашёл) | ближайший прецедент — Bremner 2001 (ранг 0 убивает 1-парам. семейство) |
| 6 | Открытая реализация CTP над ℚ(√D) с полным 2-кручением | **новое инструментально** (упоминаний не нашёл) | PARI 2.17.4 — только ℚ (проверено); Sage — нет (проверено); Magma/Donnelly — закрытая |

**Сводно:** из шести пунктов проверку на новизну уверенно выдерживает **один** (п. 5) плюс **один инструментальный**
(п. 6). Три пункта (1, 2 в основной части, 3) — переоткрытия; п. 4 — известный метод в новом частном применении.
Плюс есть «серая зона» (2d, 3′, 4″), где я упоминаний не нашёл, но это мелкие технические следствия, а не результаты.

---

## Найденные попутно вещи, которые стоит проверить отдельно

1. **Семейство A. Wesołowski (20.12.2019, multimagie.com) совпадает по форме с нашим G1.**
   Его шесть автоматических квадратов и три свободные клетки — x²+y²z², (x²+y²)(z²+1)/2, x²z²+y² — это буквально
   наши F₀, F₄, F₈ из `INTEGRALS_BRIDGE_FOR_CLAUDE_2026-09-12.md`, §3. Надо выяснить, одно и то же это семейство
   или нет; если да, то «основная новая конструкция» §2.2 в `BREMNER16_2026-09-12.md` — тоже переоткрытие.
   Бойер там же пишет: «When n is an even positive number, **we can prove that the cells in red can't be squared
   integers**» — то есть у них уже есть частичный результат о невозможности внутри этого семейства.
2. **R. Rathbun (11–13.03.2017)** нашёл эллиптическую кривую ранга 4
   y² = x³ − x² + 7528505392x + 671534074163712, дающую известный пример с 7 квадратами; «thousands of points
   on this curve … did not produce any other magic square». Это прямая перекличка с нашим MW-поиском.
3. **Ауэль–Сингер не делают спуска.** Грепы по исходнику TeX: слов «Chabauty», «Selmer», «descent»,
   «Mordell–Weil rank» (кроме ссылки на Faltings) нет. Значит арифметическая (а не геометрическая) часть
   у них действительно не закрыта — это и есть наша ниша.
4. **Формула Бремнера для степеней 16 конфигураций** (12, 8, 12, 20, 12, 12, 12, 20, 20, 12, 20, 20, 8, 12, 20, 12)
   даёт роды кривых седьмой клетки напрямую (⌈d/2⌉−1 = 5, 3, 9). Стоит сверить с нашей таблицей в `BREMNER16`.

---

## Черновик письма (английский; без завышений)

> Dear Professor Auel,
>
> I have been working on the 3×3 magic square of squares from the arithmetic side, and your paper with Benjamin Singer
> (arXiv:2609.09351) is the closest thing I have found to the geometric picture I needed, so I hope a short question
> is not an imposition. My line of attack is to fix the ratio of one complementary pair of cells — equivalently a
> primitive arithmetic progression b², n², h² with b²+h² = 2n² — and to rule that ratio out at *every* scale, rather
> than to search centres up to a bound as Buell, Boyer and Morgenstern have done. Each such ratio gives a genus-2
> curve whose Jacobian is a Weil restriction Res_{k/ℚ}E₁ for a real quadratic k, and for three ratios, (17:7:13),
> (23:7:17) and (71:49:61), I can prove rank E₁(k) = 1 and then close the curve by elliptic Chabauty; the rank-1 proof
> for k = ℚ(√165) needed the Cassels–Tate pairing on Sel²(E₁/k), which I had to implement from scratch in Sage since
> I could find no open implementation over a number field. I am aware that most of the surrounding structure is not
> new — the two-complementary-pairs reduction is already in Boyer's 2004 search note and in Morgenstern's 2010 method,
> the base conic and the link to y² = x³ − x go back to Robertson and to Bremner's 1999 paper, and the genera 0, 1, 5,
> 17, 49 that I computed for the covering tower are exactly the ones in your Propositions 3.1, 3.6 and 3.7 — so my
> question is narrow: is the ratio-by-ratio approach something you would expect to say anything, and do the genus-5
> curves I get (double covers of a genus-1 base, three quadratic conditions u_i² = F_i(t)) coincide with your degree-8
> genus-5 components, or are they genuinely different curves on V?
>
> I would be glad to send the computations and the code; I would equally value being told that the approach is a dead
> end, and why.

(Восемь предложений; при желании можно урезать до пяти, выбросив перечисление «что не ново» и вопрос про совпадение
кривых, но именно эти два места делают письмо честным.)

---

## Чего найти не удалось (сводно)

- Полного текста D. Buell, *A search for a magic hourglass* — PDF с multimagie.com не распознаётся `pdftotext`
  (сканированный/битый шрифт); содержание известно только по пересказам Бремнера, Бойера и Пьерра–Тирье–Циммермана.
- Страницы T. Fisher `dpmms.cam.ac.uk/~taf1000/progs.html` — 403 из этой среды; состав его Magma-пакетов не проверен.
- Полного текста L. Rabern, Rose-Hulman Undergrad. Math. J. 4 (2003) — сервер отдаёт 403; содержание известно
  по пересказу Моргенштерна (июль 2015), который заявляет, что переdoказывает все его свойства элементарно.
- Статьи C. Boyer, Math. Intelligencer 27 (2005), 52–64 целиком — только начало на сайте и supplement (4 стр.);
  сама статья за paywall. Ключевой для нас фрагмент (классификация 7.I–7.VIII и «two, three or four lines through
  the centre») есть в открытом Search.pdf 2004 г., который в статью и вошёл как ссылка [8].
- Открытых репозиториев с реализацией спаривания Касселса–Тейта (любого) — ни одного.
- Публикации, где бы к магическому квадрату квадратов применялись 2-спуск, группы Селмера, elliptic curve Chabauty
  или квадратичный Шаботи. Похоже, таких нет вовсе — ни у Бремнера (он работает с рангами над ℚ(λ) и с решётками
  Нерона–Севери K3), ни у Ауэля–Сингера.
