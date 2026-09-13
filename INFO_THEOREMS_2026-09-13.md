СТАТУС: литература — точные цитаты и «не найдено»; сам ничего не считал
ДЛЯ: всем; Claude и Codex — п. 1 и п. 2 (от них зависит, есть ли смысл искать бесконечные серии наклонов); Fable — п. 3–4
ИТОГ: (1) у Bruin–Thomas–Várilly-Alvarado доказана только КОНЕЧНОСТЬ числа кривых рода 0 и 1 на всей 9-клеточной поверхности; явного списка нет, авторы прямо пишут, что для неё это вне вычислительных возможностей. (2) Критерия типа Монского для y² = x(x+m²)(x+n²) не нашёл; все теоремы о распределении 2-Селмера (Swinnerton-Dyer, Kane, Smith) требуют ОТСУТСТВИЯ рациональной циклической подгруппы порядка 4, а у наших кривых кручение ℤ/2×ℤ/4. (3) Безусловного результата «бесконечно много простых ℓ, у которых нечётные простые делители ℓ²−1 все ≡ 3 (mod 4)» не нашёл; ближайший аналог с двумя сдвигами (Friedlander–Iwaniec) доказан лишь при гипотезе Эллиотта–Халберстама. (4) Разреза по наклону q = kp в литературе не нашёл, списков исключённых наклонов нет. (5) Для односвязных поверхностей общего типа с q = 0 безусловной незарисской-плотности нет; Lawrence–Venkatesh неприменим.
ОТМЕНЯЕТ: «p_g = 1» в контексте задания. У Auel–Singer, Theorem 2: h^{2,0} = 111, то есть p_g = 111. Также отменяет надежду получить из BTVA явный список кривых рода ≤ 1.
ПРОВЕРЕНО: полные тексты (pdftotext): BTVA arXiv:1912.08908v3; Bremner 1999 и 2001 (сканы с multimagie.com); Kane arXiv:1009.1365; Nath–Xie arXiv:2501.16723v3; Green–Tao arXiv:math/0606088v2 (§1); Stoll–Testa arXiv:1009.0388v2 (конец §); Rome–Yamagishi (открытый PDF, введение); Xie arXiv:2607.02165 (§1.3); Ghale–Islam arXiv:2308.00625v2 (введение); Yoshida arXiv:2407.09825 (введение). Аннотации: Smith 1702.02325, Naskręcki 1210.6933, Chan–Hanselman–Li 1805.10709, Bhargava–Ho 2207.03309, Lawrence–Venkatesh 1807.02721, Chakraborty–Ghale 2301.08099, FIMR 1110.6331 (только поисковая выдача). Auel–Singer 2609.09351 — HTML v1 через пересказ WebFetch плюс цитаты из LITSEARCH_2026-09-12.md §3.
ЧИТАЛ: LITSEARCH_2026-09-12.md §3, NOVELTY_CHECK_2026-09-12.md §5, начало REPORT_LINEAR_SECTIONS_FEASIBILITY_2026-09-13.md
ЧИТАТЬ ДАЛЬШЕ, ЕСЛИ: думаете о бесконечной серии закрытых наклонов, о «списке исключительных кривых» или хотите сослаться на литературу в письме

# Теоремы из литературы для линейных сечений (13.09.2026)

Метки: **[ПРОЧИТАЛ]** — прочитал в полном тексте, цитата дословная; **[АННОТАЦИЯ]** — видел только аннотацию или страницу arXiv;
**[ВТОРИЧНО]** — знаю только по цитированию в другой статье; **[НЕ НАЙДЕНО]**; **[НАБЛЮДЕНИЕ]** — мой вывод, в источнике его нет.

---

## 1. Bruin–Thomas–Várilly-Alvarado: что именно доказано

**Источник.** N. Bruin, J. Thomas, A. Várilly-Alvarado, «Explicit computation of symmetric differentials and its application to
quasi-hyperbolicity», *Algebra & Number Theory* **16** (2022), no. 6, 1377–1405, DOI 10.2140/ant.2022.16.1377 (выходные данные сверены на msp.org);
arXiv:1912.08908 (v3 от 12.10.2021). **[ПРОЧИТАЛ]** введение, §6 и §7.

**Определение** (§1, первый абзац): «In this article we call a complex projective surface Y *algebraically quasi-hyperbolic* if it contains
only finitely many curves of geometric genus 0 or 1.» Речь о комплексной поверхности: считаются все кривые над ℂ, а не только над ℚ.
**Какая поверхность** (§1.4.2): X_ms ⊂ P⁸ задана равенствами «x₁²+x₂²+x₃² = x₄²+x₅²+x₆² = x₇²+x₈²+x₉² = x₁²+x₄²+x₇² = x₂²+x₅²+x₈² =
x₁²+x₅²+x₉² = x₃²+x₅²+x₇²». Это **все девять клеток** (равенство третьего столбца следует из остальных), та же поверхность V, что у Auel–Singer.
«This surface is smooth except for 256 isolated ordinary double points. This exceeds ℓ_min(8) = 217, so we obtain the following result.
**Theorem 1.5.** The complex projective surface X_ms ⊂ P⁸ that parametrizes 3×3 magic squares of squares is algebraically quasi-hyperbolic.»
Метод: оценка снизу для h⁰(Y, SᵐΩ¹_Y) (Theorem 1.4) дает симметрические дифференциалы. Затем Proposition 6.2: «there are only finitely many
complete curves C on Y of genus at most 1 on which ω restricts to 0» (доказательство через теорему Жуанолу о слоениях).
Числа из статьи: «for m ≥ 47 there are global sections and that H⁰(Y_ms, S⁴⁷Ω¹) ≥ 8448».

**Явного исключительного множества нет.** §7, сразу после Theorem 1.5, цитата: «Unfortunately, X_ms is out of range of current computational
techniques to explicitly determine Ŝᵐ Ω¹_{X_ms}, so we cannot apply the methods from Corollary 6.5 get an explicit description of the locus of
special curves.» Для кубоида (Thm 1.2) и секстики Барта (Thm 1.1) явные ограничения получены: например, кривая рода 0, отличная от 32 коник
ван Луйка, проходит не менее чем через 7 особых точек. **Для магической поверхности ничего подобного нет: ни степени, ни числа особых точек,
ни списка.** Auel–Singer (arXiv:2609.09351, Remark A.2, цитата из LITSEARCH §3.2) тоже оставляют открытым даже вопрос о кониках: «We leave open
the question of whether these are all the conics contained in V; we know there are finitely many by [12]». **Ответ на вопрос 1:** утверждение
«все кривые рода ≤ 1 на V лежат в явном конечном списке и все вырождены» из литературы **не следует**. Доказана только конечность.
**[НАБЛЮДЕНИЕ]** Наши сечения при фиксированном k — кривые большого рода, а E_{m,n} — их факторы. Теорема 1.5 про такие кривые ничего не говорит.

> Мелкое расхождение. Rome–Yamagishi (Res. Number Theory 11 (2025), art. 91, введение) пишут о «lines parametrising repeated entries» на V.
> Но у Auel–Singer, Prop. A.1: «The Fano scheme of lines F(V) … is empty». Верить Auel–Singer: это теорема с доказательством, а у Rome–Yamagishi — попутная фраза.

---

## 2. Семья y² = x(x+m²)(x+n²): ранг 0, 2-Селмер, критерии

**Главное препятствие [ПРОЧИТАЛ].** Все известные теоремы о распределении ранга 2-Селмера с критериями через символы Лежандра
между простыми делителями (Монский/Heath-Brown для конгруэнтных чисел, их обобщения) относятся к **квадратичным твистам фиксированной
кривой** и требуют отсутствия циклической 4-подгруппы:
- Kane, *Algebra Number Theory* 7 (2013) 1253–1279, arXiv:1009.1365, с. 1: «E is an elliptic curve over Q with complete 2-torsion and
  **no cyclic subgroup of order 4 defined over Q**». Там же, §1: «if someone were to obtain a Swinnerton-Dyer type result for twists of an
  elliptic curve with full 2-torsion over Q that has a rational 4-isogeny, it is almost certain that the techniques from this paper would allow…»
  То есть случай с 4-изогенией тогда оставался открытым.
- Smith, arXiv:1702.02325, аннотация **[АННОТАЦИЯ]**: «Given an elliptic curve E/Q with full rational 2-torsion and **no rational cyclic
  subgroup of order four**, we analogously prove that the 2^∞-Selmer groups of the quadratic twists of E have distribution as predicted by Delaunay's heuristic.»

У E_{m,n} есть точка порядка 4, значит, она **вне гипотез** обеих работ. Кроме того, наша семья — не твист-семья: j меняется вместе с m/n.
**[НАБЛЮДЕНИЕ]** Деля на m², получаем y² = x(x+1)(x+t²), t = n/m. Это однопараметрическое семейство с неизотривиальным j, и теоремы о твистах к нему не относятся.

**Родственные семьи [АННОТАЦИЯ / ПРОЧИТАЛ введение]:**
- Ono, «Euler's concordant forms», *Acta Arith.* 78 (1996) 101–123: E_{M,N}: y² = x(x+M)(x+N). У нас M = m², N = n². Полный текст получить не удалось
  (eudml 403, matwbn не отдал PDF). Что там есть про торсию и ранг при квадратных M, N — **не проверено**.
- Naskręcki, *Acta Arith.* 160 (2013) 159–183, arXiv:1210.6933 [АННОТАЦИЯ]: «the family y² = x(x−a²)(x−b²) parametrized by Pythagorean
  triples (a,b,c). We prove that for a generic triple the lower bound of the rank … is 1». Это оценки ранга **снизу** (≥ 1), противоположное нужному нам.
- Heron/θ-конгруэнтные: Ghale–Islam arXiv:2308.00625v2 [ПРОЧИТАЛ введение] и Chakraborty–Ghale arXiv:2301.08099 [АННОТАЦИЯ].
  Там явные формулы 2-Селмера для E: y² = x(x−1)(x+n²) при условии «n²+1 = 2q, q простое», с дихотомией по классам простых делителей n mod 8.
  Форма похожа на нашу (одно из смещений — квадрат), но второе смещение равно 1, а не квадрату, и условие на q узкое. Прямо не переносится.
- Body cuboid: связь «y² = x(x+a²)(x+b²), нужна точка вне кручения» я видел **только в сниппете поисковика**. В лекции Knill
  (abel.math.harvard.edu/~knill/various/eulercuboid/lecture.pdf) этой формулы нет — **не подтверждено**. В любом случае интерес там к рангу > 0. Yoshida arXiv:2407.09825 доказывает бесконечность s с рангом > 0 для y² = x(x−(2s)²)(x+(s²−1)²).
- Chan–Hanselman–Li arXiv:1805.10709 [АННОТАЦИЯ]: семья ℤ/2×ℤ/8 (подсемейство ℤ/2×ℤ/4). Ранги получены условно (GRH+BSD), средние размеры 2-Селмера и произведения Тамагавы стремятся к ∞.
- Bhargava–Ho arXiv:2207.03309 [АННОТАЦИЯ]: средние 2- и 3-Селмера «in various families of elliptic curves with marked points».
  Есть ли среди них ℤ/2×ℤ/4 — по аннотации не видно; **не проверено**.

**[НЕ НАЙДЕНО]:** критерий типа Монского для y² = x(x+m²)(x+n²) через простые делители m, n, n²−m²; бесконечная серия этой формы с
**доказанным** рангом 0; теорема о положительной доле ранга 0 в семье ℤ/2×ℤ/4.

---

## 3. Бесконечные серии через простые

**Условие «нечётные простые делители ℓ²−1 все ≡ 3 (mod 4)»: [НЕ НАЙДЕНО].** Ближайшие безусловные и условные результаты
(все [ПРОЧИТАЛ] во введении Nath–Xie, arXiv:2501.16723v3, DOI 10.4064/aa250227-19-11, *Acta Arith.*):
- Один сдвиг, полуразмерное решето. Iwaniec, *Acta Arith.* 21 (1972) 203–234: «#{p ≤ x : p = m²+n²+1} ≫ x/(log x)^{3/2}». До этого Linnik (1960) доказал бесконечность. **[ВТОРИЧНО]**
- Два сдвига одновременно. Friedlander–Iwaniec, «Hyperbolic prime number theorem», *Acta Math.* 202 (2009) 1–19: «Assuming the
  Elliott-Halberstam Conjecture, they showed that x ≪ Σ_{n≤x} Λ(n) r(n−2) r(n+2) ≪ x» (сноска: «The upper bound is established without any assumption»). **[ВТОРИЧНО]**
- Один сдвиг точно, второй — почти простое. Nath–Xie, Theorem 1.1: «#{p ≤ x : p = m²+n²+1, Ω(p+2) ≤ 9} ≫ x/(log x)^{5/2}».

**[НАБЛЮДЕНИЕ]** Наше условие — два полуразмерных условия на соседние сдвиги p−1 и p+1 (исключены простые ≡ 1 mod 4). По сложности это
та же задача, что у Friedlander–Iwaniec, где «сумма двух квадратов» ставится на p−2 и p+2. Там даже порядок величины получен только при EH.
Ждать безусловного результата для ℓ²−1 не приходится.

**Грин–Тао–Циглер [ПРОЧИТАЛ §1 Green–Tao].** Green–Tao, «Linear equations in primes», arXiv:math/0606088 (*Ann. of Math.* 171 (2010);
безусловность при любой конечной сложности — после GTZ, arXiv:1009.3998, и работы о нильпоследовательностях). Lemma 1.6: «this system has finite complexity
if and only if no two of the ψᵢ are affinely dependent». **[НАБЛЮДЕНИЕ]** Если r, s — нечётные простые, то s±r чётны. Требовать простоты
всех четырёх форм r, s, s−r, s+r нельзя из-за локального препятствия при p = 2; остаётся r = 2, и тогда s−2, s, s+2 аффинно зависимы
(бесконечная сложность, тип близнецов), так что теорема не применяется. Для систем, где простыми требуются лишь некоторые формы, ГТЦ даёт асимптотику,
но **условий на символы Лежандра между значениями форм в этих работах нет**. Работ «линейные формы в простых + предписанные символы
Лежандра между ними» **[НЕ НАЙДЕНО]**. Ближе всего: Friedlander–Iwaniec–Mazur–Rubin, «The spin of prime ideals», *Invent. Math.* 193 (2013) 697–749,
arXiv:1110.6331 (билинейные суммы с символами над простыми, приложение к Селмеру) **[видел только выдачу поиска]**; комбинаторика
«Legendre specifications» в разборе Smith у Koymans (слайды «Smith explained III») **[видел только сниппет]**. Обе — про мультипликативные
произведения простых, не про линейные формы. Heath-Brown, *Invent. Math.* 111 (1993) 171–195 и 118 (1994) 331–370 — моменты 2-Селмера
для y² = x³ − D²x **[ВТОРИЧНО, по введению Ghale–Islam]**. Относится к твистам кривой без 4-кручения, см. п. 2.

---

## 4. Магические квадраты из квадратов через разрезы

**Параметризация есть, разреза по наклону нет.** Bremner, «On squares of squares», *Acta Arith.* 88 (1999) 289–297, **[ПРОЧИТАЛ]** §1:
«Any three-by-three magic square of rational numbers has the form [a−b, a+b+c, a−c; a+b−c, a, a−b+c; a+c, a−b−c, a+b]». Это ровно
наш нормированный вид: a = 1, p = b, q = c. «The square is trivial … precisely when bc(b²−c²)(b²−4c²)(4b²−c²) = 0». Бремнер фиксирует **c** (класс
твиста), а не отношение c/b: «Associate … the elliptic curve E: y² = x(x² − c²). … the existence of a magic square of squares is
equivalent to the existence of three points in 2E(Q) with x-coordinates in arithmetic progression. This observation appears first to have
been noticed by Robertson [6]» (Robertson, *Math. Mag.* 69 (1996) 289–293 — **[ВТОРИЧНО]**, полный текст не открывал).
Дальше у него эллиптическая поверхность над ℚ(λ) и параметрические «squared squares» с одной плохой диагональю, а не исключение параметров.

Bremner, «On squares of squares II», *Acta Arith.* 99 (2001) 289–308, **[ПРОЧИТАЛ]** §1–2: 16 конфигураций шести квадратов, эллиптические
расслоения над ℚ(λ). Единственное исключение семейства через ранг 0 (Remark после (10)): «this elliptic curve of conductor 48 has rank 0,
forcing μ = 0 and a consequent trivial square». Это семья седьмой клетки в конфигурации II, а не наклон. См. также NOVELTY_CHECK §4–5.
Rome–Yamagishi (круговой метод, n ≥ 4), Auel–Singer (геометрия V), Wolird arXiv:2310.12164, Cain arXiv:1908.03236 [АННОТАЦИИ] разрезов не делают.
«Kim», «Garcia» в связи с 3×3 из квадратов **[НЕ НАЙДЕНО]** (один поиск; возможно, путаница с García-Fritz–Urzúa, которые работают с кубоидом).

**[НЕ НАЙДЕНО]:** срез q = kp; 12 троек/четвёрок клеток с кривыми ℤ/2×ℤ/4; какие-либо списки исключённых наклонов.

---

## 5. Незарисская плотность для поверхностей общего типа с q = 0

**Безусловных теорем нет [ПРОЧИТАЛ].** Xie, arXiv:2607.02165, §1.3: «Over number fields, the foundational results are Faltings' theorem for
curves [Fal83] and Faltings' proof of the Mordell–Lang theorem for subvarieties of abelian varieties [Fal94]. … Beyond these cases, the general
non-density statement for varieties of general type … is still largely open.» Stoll–Testa (arXiv:1009.0388v2, конец статьи) о поверхности кубоида —
тоже полном пересечении квадрик с узлами. Там про **конечность кривых рода ≤ 1** сказано: «this is only known in very few cases, for instance, when Bogomolov's
inequality (K_S)² > c₂(S) holds; for surfaces contained in an abelian variety; for very general surfaces of large degree in projective space…
the surface S is simply connected and hence is not contained in an abelian variety». Для V (Auel–Singer Thm 2: K² = 576 < 768 = c₂,
односвязна, q = 0) не работает ни Богомолов, ни Фальтингс. Конечность кривых для V даёт только BTVA (п. 1). Для **точек** ничего нет.

**Lawrence–Venkatesh [АННОТАЦИЯ + НАБЛЮДЕНИЕ].** arXiv:1807.02721: новое доказательство Морделла и незарисская плотность гиперповерхностей
«with good reduction away from a fixed set of primes» **в пространстве модулей**. Метод требует семейства над базой с большой монодромией
и квазиконечным отображением периодов. Над всей Ṽ монодромия тривиальна (π₁ = 1), а подходящее семейство над открытым U никому не известно.
Поэтому метод к рациональным точкам V неприменим.

---

## Не делать

- Не ссылаться на BTVA как на «явный список исключительных кривых» или «все кривые рода ≤ 1 вырождены»: доказана только конечность (§7, цитата выше).
- Не писать p_g = 1: у Auel–Singer p_g = h^{2,0} = 111.
- Не ссылаться на Smith / Kane / Swinnerton-Dyer / Heath-Brown для E_{m,n}: гипотеза «нет циклической 4-подгруппы» нарушена, и семья не твист-семья.
- Не искать дальше «Kim» и «Garcia» как авторов исключения наклонов: следов нет.
- Не рассчитывать на Грина–Тао для четвёрки r, s, s−r, s+r, где все четыре простые: мешает чётность, остаётся система типа близнецов.
- Не обещать безусловную бесконечную серию простых с условием на ℓ²−1: даже аналог с суммами двух квадратов доказан только при EH.

## Недоступно / не проверено

- Полный текст Ono 1996 (торсия и ранг E_{M,N} при квадратных M, N) — не получил.
- Robertson 1996 — только через цитату у Бремнера.
- Полный текст Auel–Singer в этой сессии — через пересказ WebFetch; цитаты Remark A.2 и Prop. A.1 взяты из LITSEARCH_2026-09-12.md §3.
- Bhargava–Ho: есть ли семья ℤ/2×ℤ/4 — не смотрел внутрь.
