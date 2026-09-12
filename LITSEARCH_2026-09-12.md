# Литературный поиск, 2026-09-12 (Claude)

Задание: (1) первоисточник формулы Касселса для CT-спаривания при полном 2-кручении и сверка с нашей;
(2) квадратичный Шаботи для рода 2 с вещественным умножением; (3) статья Auel–Singer о поверхности магического квадрата.

Статусы: **цитата** — дословно из источника; **вывод** — моя выкладка на основе цитат; **не найдено** — честно не найдено.

---

## 1. Формула Касселса (Cassels 1998) и сверка с нашей

### 1.1. Первоисточник и доступность

**J. W. S. Cassels, «Second descents for elliptic curves», Journal für die reine und angewandte Mathematik 494 (1998), 101–127.**
(Посвящена Мартину Кнезеру к 70-летию.)

**Полного текста получить не удалось.** Том 494 у de Gruyter закрыт подпиской; DigiZeitschriften прекратил работу 31.12.2025
(открытый контент обещан перенести в репозиторий SUB Göttingen «к середине 2026»), GDZ том 494 не отдаёт; открытого скана/препринта нет
(arXiv 1998 г. по этой теме отсутствует). Поэтому ниже — **два независимых пересказа**, которые совпадают друг с другом
и явно ссылаются на §7 Касселса.

Номера у Касселса, извлекаемые из цитирований:
- **[Cas98, Lemma 7.4]** — определение спаривания (ссылка Shukla–Stoll);
- **[Cas98, Lemma 7.2]** — локальное зануление (ссылка Zhang);
- **[Cas98, §2]** — модели 2-накрытий.

### 1.2. Пересказ A: Shukla–Stoll (над произвольным числовым полем K)

H. Shukla, M. Stoll, «The Cassels–Tate pairing on 2-Selmer groups of elliptic curves», arXiv:2302.01640 (3 Feb 2023), §3.2.

Установка (цитата, (3.9)):
> «Let a, a′ ∈ S⁽²⁾(E) be represented by β = (β₁,β₂,β₃) and β′ = (β′₁,β′₂,β′₃) … define quadratic forms in variables (U₁,U₂,U₃,T) by
> H_i(U₁,U₂,U₃,T) := (β_jΓ_j² − β_kΓ_k²)/(e_j − e_k) + T², where Γ_i := U₁ + U₂e_i + U₃e_i². … Any two of these quadratic forms
> define a projective curve D_a in P³ … D_a has points for every completion K_v of K and is a 2-covering of E representing a.
> For details see [Cas98, §2].»

Точки на кониках и касательные (цитата):
> «This implies that for each i there is a point q_i := (u₁ : u₂ : u₃ : 1) or (Γ*_j : Γ*_k : 1) defined over K(e_i) satisfying H_i = 0,
> which is a consequence of the local-global principle for quadratic forms. Let L_i(U₁,U₂,U₃,T) for 1 ≤ i ≤ 3 be a linear form such that
> L_i = 0 is the tangent to H_i at q_i.»

**Формула Касселса (3.10), цитата:**
> «If q_v is a point defined on D_a over K_v, then Cassels' pairing [Cas98, Lemma 7.4] is defined as follows:
> ⟨α, α′⟩_Cas := ∏_v ∏°_i (L_i(q_v), β′_i)_{K_v(e_i)}.
> In [Cas98] Cassels showed that the above definition gives a well-defined pairing and is independent of the choices made.»

(∏°_i — произведение по одному представителю каждой G_K-орбиты на Δ = E[2]∖{T₀}; при полном K-рациональном 2-кручении это i = 1,2,3,
и K_v(e_i) = K_v.)

Равенство с CT-спариванием — **Theorem 5.6** (цитата): «For a, a′ ∈ S⁽²⁾(E), we have ⟨a,a′⟩_Cas = ⟨a,a′⟩_CT.»
Раньше это доказано в T. Fisher, E. F. Schaefer, M. Stoll, «The yoga of the Cassels–Tate pairing», LMS J. Comput. Math. 13 (2010), 451–460
(arXiv:0710.2079), аннотация: «Cassels has described a pairing on the 2-Selmer group of an elliptic curve which shares some properties
with the Cassels-Tate pairing. In this article, we prove that the two pairings are the same.»

Полезное для нас (обрезание произведения) — **Corollary 5.5**: вклад дают лишь конечные места из
S_{a,a′} ⊂ {простые плохой редукции E/K} ∪ {2}, где локализация коцикла α или α′ проходит через разветвлённое расширение;
бесконечные места вклада не дают.

### 1.3. Пересказ B: Shenxing Zhang (над ℚ, полное 2-кручение) — самый явный

S. Zhang, «On a comparison of Cassels pairings of different elliptic curves», arXiv:2303.05059, §2.1.

E: y² = x(x − e₁)(x + e₂), Sel₂(E/ℚ) ∋ Λ = (d₁,d₂,d₃), d₁d₂d₃ ≡ 1 mod ℚ^{×2}, однородное пространство (цитата, (2.1)):
> «H₁ : e₁t² + d₂u₂² − d₃u₃² = 0, H₂ : e₂t² + d₃u₃² − d₁u₁² = 0, H₃ : e₃t² + d₁u₁² − d₂u₂² = 0.»

**Формула (цитата, с. 4):**
> «Cassels in [Cas98] defined a skew-symmetric bilinear pairing ⟨−,−⟩ on the F₂-vector space Sel′₂(E). … For any Λ ∈ Sel₂(E),
> choose P = (P_v)_v ∈ D_Λ(A_ℚ). Since H_i is locally solvable everywhere, there exists Q_i ∈ H_i(ℚ) by Hasse–Minkowski principle.
> Let L_i be a linear form in three variables such that L_i = 0 is the tangent plane of H_i at Q_i. For any Λ′ = (d′₁,d′₂,d′₃) ∈ Sel₂(E), define
> ⟨Λ,Λ′⟩_E = Σ_v ⟨Λ,Λ′⟩_{E,v} ∈ F₂, where ⟨Λ,Λ′⟩_{E,v} = Σ_{i=1}³ [L_i(P_p), d′_i]_v.
> This pairing is independent of the choice of P and Q_i, and is trivial on E(ℚ)[2].»

Здесь Sel′₂(E) = Sel₂(E/ℚ)/(E(ℚ)_tors/2E(ℚ)_tors) — **«чистая» 2-Селмерова группа**; спаривание тривиально на E(ℚ)[2].

**Локальное зануление — Lemma 2.2 (цитата, ссылка [Cas98, Lemma 7.2]):**
> «The local Cassels pairing ⟨Λ,Λ′⟩_{E,p} = 0 if p ∤ 2∞, the coefficients of H_i and L_i are all integral at p, and modulo D_Λ and L_i by p,
> they define a curve of genus 1 over F_p together with tangents to it.»

**Никаких поправочных множителей ни в одном из двух пересказов нет** — ни при v = 2, ни при v | ∞.

### 1.4. Сверка с нашей формулой

Наша (REPLY_TO_ASTRA, п. 4; `descent/ctp.sage`):

  ⟨ε,η⟩ = ∏_v (f_{T₁}(P_v), η₂)_v · (f_{T₂}(P_v), η₁)_v,   f_{T_i} = L_i / L₀,

где L₀ — касательная плоскость четвёртого конуса пучка (того, в котором нет координаты u₀), L_i — конуса без u_i.

Касселс (обозначая η_i = β′_i = d′_i):

  ⟨ε,η⟩ = ∏_v (L₁(P_v), η₁)_v (L₂(P_v), η₂)_v (L₃(P_v), η₃)_v.

Совпадение обозначений конусов: у Касселса H_i содержит Γ_j, Γ_k, T, то есть **не содержит u_i** — ровно наш Q_i.
Четвёртый конус (без u₀ = T) Касселс **не использует**, у нас он входит как нормировка L₀.
Ещё различие — **перестановка индексов** (мы спариваем L₁ с η₂, L₂ с η₁; Касселс — L_i с η_i).

**Вывод (моя выкладка, статус: проверено алгебраически, не цитата):** формулы **тождественно равны**. Именно:

- Обозначим C = ∏_i (L_i,η_i), O = (L₁/L₀,η₂)(L₂/L₀,η₁). Тогда покомпонентно
  C/O = (L₁,η₁η₂)(L₂,η₁η₂)(L₃,η₃)(L₀,η₁η₂) = (L₀L₁L₂L₃(P_v), η₃)_v,
  так как η₁η₂ ≡ η₃ mod квадратов (условие Селмера η₁η₂η₃ ∈ K*²).
- Положим Λ := L₀L₁L₂L₃. На 2-накрытии D = C_ε имеем div(L_j) = 2D_j, где D₀ = D_O, D_i = D_{T_i} (Лемма 3.1 Yan).
  Отображение T ↦ [D_T − D_O] — гомоморфизм E[2] → Pic⁰(D), поэтому [D₁ + D₂ + D₃ − 3D₀] = [D_{T₁+T₂+T₃} − D_O] = 0,
  то есть D₁+D₂+D₃−3D₀ = div(h) для некоторой h ∈ K(D)* (K-рациональность — по Гильберту 90, дивизор K-рационален).
  Значит Λ / L₀⁴ = c·h² с **глобальной константой** c ∈ K*.
- Отсюда (Λ(P_v), η₃)_v = (c, η₃)_v для каждого v, и ∏_v (c,η₃)_v = 1 по формуле произведения. ⇒ ∏_v C = ∏_v O. ∎

Проверка существенности перестановки: «несвопнутая» версия N = (L₁/L₀,η₁)(L₂/L₀,η₂) даёт C/N = (f_{T₃}(P_v), η₃)_v,
а ∏_v (f_{T₃}(P_v),η₃)_v **не обязано** быть 1. То есть перестановка индексов (из спаривания Вейля) обязательна,
и наша формула её содержит правильно.

### 1.4a. Третье независимое подтверждение: общая форма Фишера–Ньютона

R. Visse, «Local computations on the Cassels–Tate pairing on an elliptic curve», магистерская, Leiden
(https://math.leidenuniv.nl/scripties/MasterVisse.pdf), Theorem 3.6 (= [FN14, Theorem 1.3]), для **нечётного** p:

> «Fact 3.4. For each T ∈ E[p](K̄), there is a degree 0 divisor a_T on C and rational functions f_T ∈ k̄(C) with sum(a_T) = T and
> div(f_T) = p·a_T. Furthermore, these f_T's may be scaled in such a way that the map (T ↦ f_T) is Galois equivariant. …
> In particular we can take a_O = 0 and f_O = 1.»
> «Theorem 3.6. … for each place v of k choose a point P_v ∈ C(k_v), avoiding the zeroes and poles of the rational functions f_T.
> Then the Cassels–Tate pairing is given by ⟨x,y⟩ = Σ_{v ∈ M_k} [f(P_v), w₁(η)]_v, which is independent of choices of P_v.»

Здесь R(C) = Map_k(E[p](k̄), k̄(C)), f = (f_T)_T, и w₁: H¹(k,E[p]) → R^×/(R^×)^p индуцировано спариванием Вейля,
т. е. **компонента w₁(η) при T есть σ ↦ e_p(η(σ), T)**. При p = 2 и полном 2-кручении это ровно η_i = x − e_i при T = T_i.

То есть «каноническая» форма — **сумма по ВСЕМ T ∈ E[2] и спаривание f_T с η_T без перестановки**:
F = ∏_v (L₁/L₀, η₁)_v (L₂/L₀, η₂)_v (L₃/L₀, η₃)_v   (f_O = 1).
Она локально (почленно!) совпадает с Касселсом: F/C = (L₀, η₁η₂η₃)_v = 1.

Наша форма — «по базису {T₁,T₂} с коэффициентами разложения», как у Yan, и именно поэтому в ней перестановка:
коэффициент n₁ при T₁ равен e₂(η(·),T₂) = η₂. По той же выкладке F/O = (c, η₃)_v с **глобальной** константой c,
поэтому ∏_v F = ∏_v O.

⚠️ **Важное для нас уточнение (статус: вывод):** наша форма и форма Касселса/Фишера–Ньютона совпадают
**только как глобальное произведение**, локальные множители различаются на (c, η₃)_v. Для `ctp.sage` это безразлично
(считается глобальное произведение по всем местам), но: (а) нельзя интерпретировать отдельный локальный символ у нас
как «локальный вклад Касселса»; (б) лемма об обрезании произведения по местам ([Cas98, Lemma 7.2]) сформулирована
для формы Касселса — наше обоснование обрезания (REPLY п. 6, и независимая проверка по 536 простым) делалось отдельно и корректно,
но формально это не цитата Касселса.

**Итог по вопросу 1: наша формула совпадает с формулой Касселса** — не «с точностью до обозначений» буквально
(у нас другое множество конусов и другая расстановка индексов), а **тождественно как глобальное произведение**;
расхождение локальных множителей равно (c, η₃)_v с глобальной константой c и гасится формулой произведения.
Дополнительно: наши наблюдения «T₁, T₂ в ядре» и «конечность множества мест» — это в точности
«trivial on E(ℚ)[2]» (Zhang, §2.1) и [Cas98, Lemma 7.2] (= наш п. 6 в REPLY_TO_ASTRA), то есть тоже согласуются с первоисточником.

### 1.5. Сверка с Yan (род 2)

J. Yan, «Computing the Cassels–Tate pairing for genus two Jacobians with rational two torsion points», arXiv:2109.08258.

**Theorem 3.3 (цитата):**
> «The Cassels-Tate pairing ⟨,⟩_CT : Sel²(J) × Sel²(J) → {±1} is given by
> ⟨ϵ,η⟩_CT = ∏_v (f_P(P_v), b)_v (f_Q(P_v), a)_v (f_R(P_v), d)_v (f_S(P_v), c)_v,
> where (,)_v denotes the Hilbert symbol for place v and P_v is an arbitrary local point on J_ϵ avoiding zeros/poles of f_P, f_Q, f_R, f_S.»

**Lemma 3.1 (цитата):** «… (i) There exists a K-rational divisor D_T on J_ϵ which represents the divisor class of φ_ϵ*(τ_{T₁}*(2Θ)).
(ii) … Suppose T is a two torsion point. Then 2D_T − 2D is a K-rational principal divisor.»

**Remark 3.2 (цитата):** «… there exist K-rational functions f_P, f_Q, f_R, f_S on J_ϵ such that div(f_T) = 2D_T − 2D for T = P,Q,R,S.»

Определение (a,b,c,d) (цитата): «Let (a,b,c,d) denote the image of η via H¹(G_K, J[2]) ≅ (K*/(K*)²)⁴ … induced by taking
the Weil pairing with P,Q,R,S», где P,Q,R,S — базис J[2].

То есть **перестановка (f_P ↔ b, f_Q ↔ a) присутствует и у Yan** — она и есть след спаривания Вейля.
Наш эллиптический аналог Thm 3.3 воспроизведён верно; и, как показано в 1.4, он эквивалентен формуле Касселса.

### 1.6. Где ещё формула Касселса пересказана (для дальнейшей сверки)

- T. Fisher, «On binary quartics and the Cassels–Tate pairing», Res. Number Theory 8 (2022), no. 74, arXiv:2208.14977.
  Цитата: «Cassels [6] also described a method for computing the pairing in the case n = 2. His method involves solving conics over
  the field of definition of each 2-torsion point on E.» и: «More recently, Donnelly [10] found a method that only involves solving conics over K.»
  Сам Фишер даёт формулу **без решения коник вообще** (инвариантная теория бинарных квартик, K3-поверхность (2,2,2)-формы) —
  это альтернатива, которую при желании можно использовать как независимую реализацию.
- S. Zhang, arXiv:2303.05059 §2.1 (см. выше) — самый явный вариант при полном 2-кручении.
- M. van Beek, T. Fisher, «Computing the Cassels–Tate pairing on 3-isogeny Selmer groups via cubic norm equations», Acta Arith. 185 (2018), 367–396.
- T. Fisher, R. Newton, «Computing the Cassels–Tate pairing on the 3-Selmer group of an elliptic curve», IJNT 10 (2014), 1881–1907.
- P. Swinnerton-Dyer, «2^n-descent on elliptic curves for all n», JLMS 87 (2013), 707–723.
- Магистерская Leiden: Visse, «Local computations on the Cassels–Tate pairing on an elliptic curve»
  (https://math.leidenuniv.nl/scripties/MasterVisse.pdf) — не разобрана, но вероятно содержит полный разбор локальных вкладов.
- H. Shukla, «Computing Cassels–Tate pairing for odd-degree hyperelliptic curves» (в подготовке, [Shu23]) — продолжение Shukla–Stoll на род ≥ 2.

### 1.7. Практическая рекомендация

Для **независимой проверки формулы** (а не только вычисления) лучший ход — реализовать вариант **Zhang §2.1 / Cassels**:
три конуса, три касательные, произведение ∏_{i=1}^{3}(L_i(P_v), η_i)_v, **без нормировки на L₀ и без перестановки**.
По п. 1.4 он обязан дать ту же матрицу Грама; если даст — это независимая перекрёстная проверка и вычисления, и самой формулы.
Это дёшево: у нас уже есть все конусы, рациональные точки на них и локальные точки —
достаточно поменять две строки в `descent/ctp.sage` (убрать деление на L₀ и снять перестановку индексов).
Заодно это снимет последнюю непроверенную зависимость от «теории Cassels/Yan»: форма Касселса — дословно цитата из двух
независимых источников, а не наш вывод.

Дополнительный бонус: в форме Касселса работает его собственная **Lemma 7.2** об обрезании произведения
(целость коэффициентов H_i и L_i в p, хорошая редукция), то есть обоснование конечного списка мест станет ссылкой,
а не самодельной леммой.

### 1.8. Готовые реализации (для независимой сверки)

- **Magma, `CasselsTatePairing`** — реализация S. Donnelly (подход через homogeneous space definition), считает спаривание
  на Sel₂(E)/E(ℚ)[2] после 2-спуска. Вопрос, работает ли она **над числовым полем** (нам нужно k = ℚ(√165)), по документации
  однозначно установить не удалось — часть источников утверждает «over number fields», проверить без Magma нельзя. **Не проверено.**
- **T. Fisher, «On binary quartics and the Cassels–Tate pairing»** (arXiv:2208.14977) — формула **без решения коник вообще**;
  как независимая реализация технически проще всего (нет уравнений нормы, нет проблемы с группой классов ℚ(√462),
  которая у нас упиралась в стек PARI).
- **T. Fisher, J. Yan, «Computing the Cassels–Tate pairing on the 2-Selmer group of a genus 2 Jacobian»**, arXiv:2306.06011
  (также https://www.dpmms.cam.ac.uk/~taf1000/papers/genus2ctp.pdf) — прямое продолжение Yan arXiv:2109.08258;
  если понадобится считать CTP на Sel²(J₂) напрямую (а не на E₁/k), смотреть сюда.
- **H. Shukla, M. Stoll**, arXiv:2302.01640 — Albanese–Albanese определение, над произвольным числовым полем;
  Theorem 5.4 + Corollary 5.5 дают вариант, в котором список мест ограничен {плохая редукция} ∪ {2} **по построению**.

---

## 2. Квадратичный Шаботи

### 2.1 (а) Точная формулировка условия и первоисточник

**J. S. Balakrishnan, N. Dogra, «Quadratic Chabauty and rational points I: p-adic heights», Duke Math. J. 167 (2018), no. 11, 1981–2038;
arXiv:1601.00388.** Неравенство r < g + ρ(J) − 1 — это **Lemma 3.2**:

> «**Lemma 3.2.** Suppose X is a curve of genus g, such that rk J(K) < g + ρ(J) − 1. Then X(K_𝔭)₂ is finite.»

Обстановка (введение), цитата:
> «Let K be **Q** or an imaginary quadratic field, and let X/K be a smooth projective curve of genus g > 1 with a K-rational point b.
> Let T₀ be the set of primes of bad reduction for X, let p be a prime of **Q** such that {v|p} ∩ T₀ is empty… Let ρ(J) = rk NS(J)
> denote the Picard number of J (**over K, not necessarily its geometric Picard number**).»

**Пять условий, которые легко потерять:**
1. ρ — ранг NS(J) **над базовым полем K**, не геометрический;
2. Lemma 3.2 даёт только **конечность**; для явного вычисления нужна **Theorem 1.2**, где гипотезы жёстче:
   «Suppose r = g, ρ(J) > 1, and the p-adic closure J(K)‾ has finite index in J(K_𝔭)»;
3. нужна рациональная точка b ∈ X(K);
4. нужна **хорошая редукция в p** ({v|p} ∩ T₀ = ∅);
5. **K = ℚ или мнимое квадратичное** (в BD I; см. таблицу в BD II).

**Balakrishnan–Dogra II**, «Quadratic Chabauty and rational points II: Generalised height functions on Selmer varieties»,
IMRN 2021, no. 15, 11923–12008; arXiv:1705.00401. Там сводная таблица: QC1 (= BD I) — «X(K) for X/K with r < g + ρ(J) − 1,
K = **Q** or im. quad.». Также Theorem 1.1 (расширение на r > g через J ~ A^d × B).

Алгоритмическая версия: **BDMTV** = Balakrishnan–Dogra–Müller–Tuitman–Vonk, «Quadratic Chabauty for modular curves:
Algorithms and examples», Compositio Math., arXiv:2101.01862, §2.1 и Algorithm 3.12:
> «Consider a smooth projective curve X_**Q** of genus g ≥ 2 whose Jacobian J has rank r = g. We also assume that the abelian
> logarithm induces an isomorphism log: J(**Q**) ⊗ **Q**_p → H⁰(X_{**Q**_p}, Ω¹)^∨ and that X(**Q**) is non-empty… Suppose that
> the Néron–Severi rank rk_**Z** NS(J) is at least 2, so that there exists a nontrivial class Z ∈ Ker(NS(J) → NS(X) ≃ **Z**).»
Плюс: «A prime p of good reduction such that the Hecke operator T_p generates End⁰(J)»; Remark 3.11: «These conditions imply that
the curve X has good reduction at p.»

**Наш случай формально проходит**: g = 2, ρ = 2, r = 2 ⇒ r = 2 < 3 = g + ρ − 1 ✓ и r = g ✓.

Альтернативный подход: **Besser–Müller–Srinivasan, «p-adic adelic metrics and Quadratic Chabauty I», Crelle 828 (2025);
arXiv:2112.03873** (через p-адическую теорию Аракелова — проще).

### 2.2 ⚠️ Главный риск: ρ(J) считается НАД ℚ

Статус: **вывод субагента, не цитата** — но проверить обязательно.

End_ℚ(Res_{k/ℚ}E₁) = Gal(k/ℚ)-инвариантная часть End_k̄(E₁ × E₁^σ).
- Если E₁ **не** изогенна сопряжённой E₁^σ над k̄ — то End⁰_ℚ(J) = ℚ, два класса E₁×0 и 0×E₁^σ в NS(J_ℚ̄) переставляются Галуа,
  и **ρ(J) = 1 над ℚ**. Тогда g + ρ − 1 = 2 и условие требует rank J(ℚ) ≤ 1 — то есть обычный Шаботи–Коулман,
  а квадратичный **не даёт ничего сверх него**.
- Если E₁ — **ℚ-кривая** (изогенна E₁^σ над k̄; это совместимо с тем, что E₁ не определена над ℚ) — тогда Res_{k/ℚ}E₁ типа GL₂,
  RM, **ρ(J) = 2**, и всё работает. Опора: J. Quer, «Q-curves and abelian varieties of GL₂-type», PLMS 81 (2000), 285–317;
  K. Ribet, «Abelian varieties over Q and modular forms», arXiv:alg-geom/9208002.
  BDMTV подтверждают связь: «The Jacobian J₀⁺(N) of X₀⁺(N) has real multiplication over **Q**, so the Picard number is at least g».

**Для нас (проверить!):** по `FAMILY_SECTIONS.md` в сечении (71,49,61) E₁ ≅ E₀ ⊗ √165 с E₀/ℚ, то есть E₁ — квадратичная закрутка
кривой над ℚ элементом √165 ∈ k. Тогда E₁^σ = E₀ ⊗ (−√165), и E₁, E₁^σ **геометрически изоморфны** (обе ≅ E₀ над k̄)
⇒ E₁ — ℚ-кривая ⇒ ρ(J) = 2 и RM, скорее всего, действительно есть. Но это надо подтвердить явно
(сравнение a_𝔭 для E₁ и E₁^σ, либо Costa–Mascot–Sijsling–Voight, «Rigorous computation of the endomorphism ring of a Jacobian»,
Math. Comp. 2019 — пакет `endomorphisms` Sijsling, Magma; частичный порт в Sage).

**Приятное следствие, если E₁ — ℚ-кривая:** тогда J типа GL₂ над ℚ ⇒ модулярна (Ribet + Khare–Wintenberger) ⇒ J — фактор J₁(N),
и работает весь модулярный аппарат QCMod, включая автоматическое получение класса Z через Эйхлера–Шимуру T_p = F_p + pF_p^{−1}.

### 2.3 (б) Реализации

**QCMod** — github.com/steffenmueller/QCMod (Balakrishnan, Dogra, Müller, Tuitman, Vonk; вклад Stoll), **Magma**.
README, цитата (условия на X/ℚ):
> «• rank = genus • **real multiplication** • enough rational points on the curve to solve for the height pairing, **unless the genus is 2**
> • p is a prime of good reduction such that – the closure of Jac(X)(**Q**) in Jac(X)(**Q**_p) has finite index
> – the Hecke operator at p generates the Hecke algebra.»
> «**QCModAffine**: Main function, takes a plane affine curve (not necessarily smooth) with integer coefficients, monic in y, and a prime p…»

**Требования к модели: не нечётная степень**, а плоская аффинная модель, моничная по y, с целыми коэффициентами,
удовлетворяющая Assumption 3.10 (Tuitman). Примеры BDMTV рода 2 — модели степени 6 (X₁₁, X₁₅), при этом для спаривания высот
они переходят к нечётностепенной модели **над ℚ_p**.

**RM-якобианы поддерживаются напрямую** (не только бисэллиптические): примеры 5.16–5.17 (X₁₁, X₁₅ — «Their Jacobians have
real multiplication, no rational torsion and Mordell–Weil rank 2; they are both absolutely simple. The Galois action on the 2-torsion
field is A₅, which is **too large for an elliptic curve Chabauty computation**») и 5.18–5.19 (C₁₈₈, C₁₆₁ — «the Jacobian of C_N is an
optimal quotient of J₀(N), so it has **real multiplication and Picard number 2**. The Mordell–Weil ranks are both 2»).

**Нужен ли явный RM-эндоморфизм?** — Remark 3.17 BDMTV, цитата:
> «*This is the only part of our algorithm specific to modular curves*, since it relies on the Eichler–Shimura relation… More generally,
> for a smooth projective curve X/**Q** satisfying the assumptions of §2.1, one could find p-adic approximations of the action of the
> nontrivial classes Z_i on H¹_dR(X_{**Q**_p}) using just p-adic linear algebra… if one knows a set of generators of a finite index
> subgroup of End(J) in advance (e.g. using algorithms for rigorous computation of the endomorphism algebra [CMSV19]) then one can use
> this to compute the classes of generators in cohomology… **When the assumption is not satisfied, our implementation throws an error,
> urging the user to try a different choice of prime p.**»
И **Remark 3.18**: «The code is currently restricted to the base field K = **Q**.»
Предупреждение авторов о немодулярных кривых: «The main stumbling block in attempting such a generalisation is our running assumption
on the Mordell–Weil rank and Picard number of the Jacobian… Since a generic curve has Picard number one…»

Локальные высоты вне p: отдельный код — **L. Betts, J. Duque-Rosero, S. Hashimoto, P. Spelier, «Local heights on hyperelliptic curves
and quadratic Chabauty», arXiv:2401.05228**, Magma: github.com/sachihashimoto/local-heights.

**Sage-ветка (единственная полноценная):** F. Bianchi, «Quadratic Chabauty for (bi)elliptic curves and Kim's conjecture»,
Algebra & Number Theory 14 (2020), 2369–2416; arXiv:**1904.04622** (не 1906.04143 — такого номера с таким содержанием нет);
код github.com/bianchifrancesca/quadratic_chabauty и github.com/bianchifrancesca/QC_bielliptic
(к Bianchi–Padurariu, «Rational points on rank 2 genus 2 bielliptic curves in the LMFDB», Contemp. Math. 796 (2024)).
README: «SageMath code for Quadratic Chabauty for rank 2, genus 2 bielliptic curves over **Q**»; кривая обязана иметь вид
y² = a₆x⁶ + a₄x⁴ + a₂x² + a₀ и «the corresponding elliptic curves each have rank 1». **Бисэллиптичность над ℚ обязательна.**

**Sage, p-адические высоты:** S. Gajović, J. S. Müller, «Computing p-adic heights on hyperelliptic curves», arXiv:2307.15787 —
«works for both odd and even degree… an implementation in SageMath». Coleman-интегрирование в Sage —
Balakrishnan–Bradshaw–Kedlaya, arXiv:1004.4936 (+ чётные степени, Balakrishnan 2015). **Полного QC-конвейера в Sage нет.**

### 2.4 (в) Бисэллиптическая ветка и Res_{k/ℚ}E

BD I **Theorem 1.4** (для y² = x⁶ + a₄x⁴ + a₂x² + a₀), цитата:
> «Suppose E₁ and E₂ each have Mordell–Weil rank 1 over K, and let P_i ∈ E_i(K) be points of infinite order… Then X(K) is contained
> in the finite set of z in X(K_𝔭) satisfying h_{E₁,𝔭}(f₁(z)) − h_{E₂,𝔭}(f₂(z)) − 2χ_𝔭(x(z)) − α₁ log_{E₁}(f₁(z))² + α₂ log_{E₂}(f₂(z))² ∈ Ω.»
BBBM-переформулировка: «Suppose that K is **Q** or an imaginary quadratic field in which p splits, and that E₁ and E₂ each have rank 1 over K…»

**Вывод: ветка требует, чтобы f₁, f₂: X → E₁, E₂ были определены над базовым полем.** Если кривая бисэллиптична только над k
(а факторы сопряжены), над ℚ этой структуры нет ⇒ **ветка неприменима**. Над k — применима, но тогда K = k вещественно квадратичное.

**QC над числовым полем:** J. S. Balakrishnan, A. Besser, F. Bianchi, J. S. Müller, «Explicit quadratic Chabauty over number fields»,
Israel J. Math. 243 (2021), 185–232; arXiv:1910.04653. Аннотация:
> «We generalize the explicit quadratic Chabauty techniques for integral points on odd degree hyperelliptic curves and for rational points
> on genus 2 bielliptic curves to **arbitrary number fields using restriction of scalars**…»

**Их бюджетное условие (1.4), цитата — прямо про нас:**
> «Over the rational numbers, the space of continuous **Q**_p-valued idele class characters has dimension 1… In general, the dimension of
> this space is at least r₂ + 1, where r₂ is the number of conjugate pairs of non-real embeddings of K into **C**… Combining our functions
> with Siksek's work, we expect that generically, we get a method to compute integral or rational points when
> (1.4)  r + rk(𝒪_K^×) ≤ [K : **Q**] · g.
> As in Siksek's work, our approach will usually fail if X can be defined over a subfield F of K and if rk(Jac(X)/F) + rk(𝒪_F^×) > [F : **Q**] · g.»

**Для вещественного квадратичного k: r₂ = 0 ⇒ всего один идельный характер (циклотомический) — над k метод даёт ровно столько же
функций, сколько над ℚ, никакого выигрыша; плюс rk(𝒪_k^×) = 1 съедает бюджет.**
Вещественно-квадратичный пример у них всё же есть — **Example 7.2**: «the **Q**(√34)-rational points on the bielliptic curve
X : y² = x⁶ + x² + 1… The rank of the Jacobian of X over **Q** is 2 and the rank over **Q**(√34) is 3.»
Реализация — **гибрид Sage + Magma**; Remark 6.6: «our current implementation is also limited to quadratic fields»;
нужны `RegularModel` (Magma) и MW-решето над K. Публичного репозитория для неё найти не удалось.

Свежее: «Rational points on the non-split Cartan modular curve of level 27 and quadratic Chabauty over number fields»,
arXiv:2501.07833 (2025) — но поле ℚ(ζ₃), мнимое квадратичное.

**Отдельной работы «quadratic Chabauty для J ~ Res_{k/ℚ}E с вещественным квадратичным k» не существует** (не найдено).

### 2.5 (г) Альтернативы при rank E(k) = [k:ℚ]

Ограничение Брюина: **N. Bruin, «Chabauty methods using elliptic curves», Crelle 562 (2003), 27–49** — метод применим только при
**rk E(K) < [K:ℚ]**. Magma handbook: «If one of the elliptic curves has rank greater than or equal to the degree of its base field,
then Chabauty's method cannot be applied» (цитата по сниппету; страница handbook отдала 401).

Что применяют вместо:
1. **Решето Морделла–Вейля** — N. Bruin, M. Stoll, «The Mordell–Weil sieve: proving non-existence of rational points on curves»,
   LMS J. Comput. Math. 13 (2010), 272–306; arXiv:0906.1934. В BDMTV применяется именно для рода 2; реализация Штолля встроена в QCMod (`mws_qc.m`).
2. **Two-cover descent** — N. Bruin, M. Stoll, «Two-cover descent on hyperelliptic curves», Math. Comp. 78 (2009), 2347–2370;
   arXiv:0803.2052 (код cecm.sfu.ca/~nbruin/twocovdesc); D. R. Hast, «Explicit two-cover descent for genus 2 curves», arXiv:2009.10313.
3. **Covering collections + elliptic curve Chabauty на накрытиях** — ровно наш тип ситуации. BDMTV §5.4, цитата:
   > «However, we note that the rational points on both curves can be computed by combining **covering collections with elliptic curve
   > Chabauty**. For C₁₈₈ this was pointed out to us by **Nils Bruin**, and for C₁₆₁, this computation is due to **Bars, González, and Xarles [BGX21]**.»
   [BGX21] = F. Bars, J. González, X. Xarles, «Hyperelliptic parametrizations of Q-curves», Ramanujan J. 56 (2021), 103–120,
   препринт **arXiv:1910.10545** — **буквально про ℚ-кривые**. Подробный разбор — в п. 2bis.7 (в частности: ранг на накрытии
   **не падает**, выигрыш даёт рост степени поля K; и наш rank E₁(k) = 2 относится к другому уровню, чем условие BGX).
   Первоисточник техники: N. Bruin, «Chabauty methods and covering techniques applied to generalized Fermat equations», CWI Tract 133 (2002).
4. **S. Siksek, «Explicit Chabauty over number fields», ANT 7 (2013), 765–793; arXiv:1010.2603** — критерий r ≤ d(g−1).
   Для g = 1 правая часть 0 (бесполезно); для рода 2 над квадратичным полем d(g−1) = 2.
5. «Extending Elliptic Curve Chabauty to Higher Genus Curves», arXiv:1111.5506 — covering techniques + модифицированное решето МВ.

Случаев «rank = степени и всё-таки elliptic curve Chabauty» не существует (это исключено устройством метода);
но есть случаи, где ранг слишком велик и задачу закрыли иначе: BDMTV Examples 5.16–5.17 (QC + решето МВ) и BGX21 (накрытия + elliptic Chabauty).

### 2.6 Что реально запустить БЕЗ Magma

| Компонент | Sage/PARI | Magma |
|---|---|---|
| Coleman-интегралы (odd & even degree) | ✅ `coleman_integrals_on_basis` | — |
| Локальная высота Коулмана–Гросса в p | ✅ Gajović–Müller, arXiv:2307.15787 | — |
| QC для бисэллиптических рода 2 **над ℚ** | ✅ `QC_bielliptic` (Bianchi) | — |
| QC для RM рода 2 общего вида (QCMod) | ❌ | ✅ `QCModAffine` |
| Подъём Фробениуса Tuitman'а | ❌ | ✅ `coleman.m` |
| Локальные высоты вне p | ❌ | ✅ `sachihashimoto/local-heights` |
| Полустабильные/регулярные модели | ⚠️ MCLF, `genus2reduction` | ✅ `RegularModel` |
| Elliptic curve Chabauty (Bruin) | ❌ | ✅ `Chabauty(...)` |
| Two-cover descent | ❌ | ✅ `TwoCoverDescent` |
| Решето МВ для рода 2 | ⚠️ частично | ✅ `mws_qc.m` (Magma ≥ 2.25) |
| Строгое End(J) (CMSV19) | ⚠️ частичный порт | ✅ `endomorphisms` |

**Честный итог: полноценный квадратичный Шаботи для рода 2 с RM без Magma сейчас запустить нельзя.**

### 2.7 Что это значит для нас

> **ОБНОВЛЕНО 2026-09-12 (см. раздел 2bis): вопрос решён отрицательно. ρ(J/ℚ) = 1, квадратичный Шаботи неприменим.**
> Пункты 1–2 ниже писались до этого расчёта и оставлены как история рассуждения; пункты 3–5 остаются в силе.

1. **Формально условие выполнено** (r = 2 < 3 = g + ρ − 1, r = g = 2) — но **только если ρ(J) = 2 над ℚ**, то есть если E₁ — ℚ-кривая.
   Это первое, что надо проверить (см. 2.2); по виду модели (закрутка E₀/ℚ элементом √D) ответ, скорее всего, положительный.
2. **Если ρ(J) = 2 подтвердится — самый перспективный путь: QCMod (Magma)** для кривой рода 2 **над ℚ**, а не для E₁ над k.
   Это ровно класс BDMTV Examples 5.16–5.19. Требования: моничная по y плоская модель с целыми коэффициентами,
   хорошая редукция в p, T_p порождает алгебру эндоморфизмов (иначе — ручная подача класса Z, Remark 3.17).
   Это задача для Codex (у него Magma).
3. **Ветка «QC над k» (BBBM) для нас почти наверняка пуста**: вещественное квадратичное поле даёт r₂ = 0 ⇒ один характер,
   плюс rk(𝒪_k^×) = 1; бюджет (1.4) не сходится. Подтверждает вывод Grok от 2026-09-11 — но по другой, более точной причине.
4. **Бисэллиптическая Sage-ветка (Bianchi) нам не подходит** — она требует бисэллиптичности над ℚ.
5. **Реалистичная альтернатива для 6 сечений с w = +1: covering collections + elliptic curve Chabauty** по образцу
   Bruin / BGX21 — и это буквально про ℚ-кривые, то есть про наш объект. Плюс решето Морделла–Вейля, которое у нас уже работает
   (`family/mw_sieve.sage`). Это, вероятно, дешевле, чем поднимать QC.
6. Ошибки в исходном задании: arXiv:1906.04143 — не та работа (нужны 1904.04622 Bianchi и 1910.04653 BBBM);
   обзор Balakrishnan «Quadratic Chabauty» в трудах ICM 2022 подтвердить не удалось.

---

## 2bis. Уточнение: степень ℚ-кривой, End и ρ(NS)

Добавлено 2026-09-12 по запросу координатора. Входные данные: у всех 11 сечений j(E₁) ∈ ℚ; для (71,49,61)
E₁ = E₀ ⊗ √165 с E₀/ℚ, так что E₁^σ = E₀ ⊗ (−√165) = E₁ ⊗ (−1), и E₁ ≅ E₁^σ **над k(i)**, но не над k.

**Короткий ответ: End⁰_ℚ(Res_{k/ℚ}E₁) = ℚ, ρ(J/ℚ) = 1, квадратичный Шаботи по Balakrishnan–Dogra
для наших сечений бесполезен. Это не «скорее всего», а следствие точной формулы Рибе.**

### 2bis.1. Формула Рибе для End⁰ ограничения Вейля

**K. A. Ribet, «Abelian varieties over Q and modular forms», Proc. KAIST Math. Workshop (1992), 53–79**
(свободный текст: math.berkeley.edu/~ribet/Articles/korea.pdf); §6, цитата:

> «Let B be the abelian variety Res_{K/**Q**} C… It represents the functor on **Q**-schemes S ↦ C(S_K); in particular, we have
> Hom_**Q**(X, B) = **Hom_K(X_K, C)** whenever X is an abelian variety over **Q**. Applying this formula in the special case where
> X = B, we find End_**Q**(B) = Hom_K(B_K, C). On the other hand, B_K = ∏_{σ∈Gal(K/**Q**)} σC. Hence we have
> **Q** ⊗ End_**Q**(B) = ∏_σ **Q** ⊗ **Hom_K(σC, C)**.»

**Ключевое слово — Hom_K: изогении между сопряжёнными должны быть рациональны над K.** Именно поэтому в §6 Рибе пишет:
> «Enlarging K if necessary, we may assume that there are isogenies μ_g : gC₀ → C₀ **defined over K**.»

(В терминологии González–Guàrdia–Rotger это «ℚ-кривая, **completely defined over K**».)

§7 «Q-curves over quadratic fields» — ровно наш случай, цитата:
> «Suppose that C is a **Q**-curve as above and that K is a quadratic field. Let σ be the non-trivial automorphism of K over **Q**.
> Then, by hypothesis, there is a **K-isogeny** μ = μ_σ : σC → C… The cocycle c takes the value 1 on all elements of
> Gal(K/**Q**) × Gal(K/**Q**) other than (σ,σ). Its value on that pair is the non-zero integer m such that μ∘σμ is multiplication by m on C.
> The algebra R may be written **Q**[X]/(X² − m)…»
> «If m is a perfect square, then we have E = **Q** in the notation of §6. The abelian variety A is then a model of C over **Q**.
> Assume for the rest of this § that m is not a perfect square. Then R = E is a quadratic number field, and we have B = A…
> **The field E is real if m is positive and imaginary if m is negative.**»

И (7.2) Proposition [Serre]: «At least one of the two quadratic fields E, K is a real quadratic field.»

Параллельная формулировка в **J. González, J. Guàrdia, V. Rotger, «Abelian surfaces of GL₂-type as Jacobians of curves»,
arXiv:math/0409352**: Definition 2.7 («A is of GL₂-type over k if End_k A is an order in a number field F of degree [F:ℚ] = dim A»)
и Lemma 3.6: если A — абелева поверхность **вещественного** GL₂-типа, k-изогенная Res_{K/k}(C), то C — ℚ-кривая,
**completely defined over K**, с j(C) ∉ k, и **F = ℚ(√(deg μ))**, где μ∘μ^σ = deg μ.

### 2bis.2. Наш случай: μ не k-рациональна ⇒ End⁰_ℚ = ℚ

E₁ не CM (j(E₁) ∈ ℚ, кондуктор E₀ = 2⁵·3²·5²·11²·13·61²·337·3061 — не CM-кривая), поэтому
Hom_k̄(E₁^σ, E₁) = ℤ·μ, где μ — изоморфизм (степень 1).

Для τ ∈ Gal(k̄/k) имеем μ^τ = ±μ (Aut(E₁) = {±1}, ибо j ≠ 0, 1728), и τ ↦ μ^τ/μ — квадратичный характер,
вырезающий k(i)/k. Он **нетривиален**, потому что −1 ∉ k*² (k = ℚ(√165) вещественно). Значит

  **Hom_k(E₁^σ, E₁) = 0**,  End_k(E₁) = ℤ,

и по формуле Рибе **ℚ ⊗ End_ℚ(Res_{k/ℚ}E₁) = ℚ ⊕ 0 = ℚ.**

Проверка через представления Галуа (независимая): V_ℓ(E₁) = V_ℓ(E₀)|_{G_k} ⊗ χ, где χ — квадратичный характер k(165^{1/4})/k;
по формуле проекции V_ℓ(Res) = V_ℓ(E₀) ⊗ Ind_{G_k}^{G_ℚ}χ. Ind χ — неприводимое двумерное диэдральное представление (D₄),
V_ℓ(E₀) неприводимо и без CM; тензорное произведение двух неприводимых двумерных представлений имеет
End_{G_ℚ} = ℚ_ℓ, и по Фалтингсу End⁰_ℚ(Res) ⊗ ℚ_ℓ = ℚ_ℓ. **Совпадает.**

Итог по вопросу 1 координатора: **не ℚ(√d) с d = 1, не ℚ(i), а просто ℚ.** Res **ℚ-прост** и **НЕ является GL₂-типом**.
Формула F = ℚ(√deg μ) здесь неприменима именно потому, что она требует k-рациональности μ; у нас μ живёт над k(i).
Разрешение парадокса «d = 1 ⇒ F = ℚ, что несовместимо с GL₂-типом размерности 2»: при d = 1 и k-рациональной μ было бы
m = ±1; при m = +1 кривая **опускается на ℚ** (Рибе: «The abelian variety A is then a model of C over **Q**»), при m = −1
получилось бы E = ℚ(i). У нас же μ вовсе не k-рациональна, и оба варианта не наступают.
Конструкция GL₂-типа для нашей ℚ-кривой требует **увеличить поле** (Рибе §6: «Enlarging K if necessary») —
над K′ = ℚ(√165, i) получается 4-мерное B = Res_{K′/ℚ}E₁ и GL₂-типный кусок A ⊆ B; но это **не наш якобиан**.

### 2bis.3. ρ(NS) над ℚ

Стандартно: при фиксированной ℚ-рациональной поляризации λ₀ (для якобиана — Θ)
NS(J) ⊗ ℚ ≅ {f ∈ End⁰(J) : f^† = f}, где † — инволюция Розати; изоморфизм Галуа-эквивариантен, поэтому

  **ρ(J/ℚ) = dim_ℚ {f ∈ End⁰_ℚ(J) : f^† = f}.**

(Mumford, *Abelian Varieties*, §20–21; Birkenhake–Lange, *Complex Abelian Varieties*, Prop. 5.2.1.)

- У нас End⁰_ℚ(J) = ℚ ⇒ симметричная часть = ℚ ⇒ **ρ(J/ℚ) = 1.**
- Контроль: End⁰_ℚ̄(J) = End⁰(E₁ × E₁) = M₂(ℚ), симметричная часть — симметричные 2×2 матрицы, размерность 3
  ⇒ ρ_geom(J) = 3 ✓ (совпадает с ρ(E × E) = 3 для неCM-кривой).
- Над k: End⁰_k(J_k) = ℚ × ℚ (диагональ, так как Hom_k(E₁,E₁^σ) = 0) ⇒ **ρ(J/k) = 2**.
- Над k(i): End⁰ = M₂(ℚ) ⇒ ρ = 3.

**Отдельное наблюдение, важное само по себе:** даже если бы получилось End⁰_ℚ = **мнимое** квадратичное поле (GL₂-тип
с CM-полем), было бы всё равно ρ = 1. Рибе, §3, дословно:
> «Let ¯ denote the canonical involution on E: this involution is the identity if E is totally real, and the “complex conjugation”
> on E if E is a CM field. **The involution is the Rosati involution on E induced by every polarization of A/Q.**»
> «Then E is a number field which is either a totally real number field or a “CM field,” since each **Q**-polarization of A defines
> a positive involution on E.»

То есть **ρ = 2 даёт только ВЕЩЕСТВЕННОЕ умножение**; при мнимом квадратичном End Розати = комплексное сопряжение,
симметричная часть = ℚ, ρ = 1. Поэтому условие Balakrishnan–Dogra и формулируется через NS, а не через End: GL₂-тип сам по себе
**не** гарантирует применимость квадратичного Шаботи — нужен именно вещественный (RM) GL₂-тип.

### 2bis.4. Вывод: квадратичный Шаботи неприменим

g + ρ(J) − 1 = 2 + 1 − 1 = **2**, а rank J(ℚ) = 2. Условие BD I, Lemma 3.2 (rank < g + ρ − 1) **нарушено**.
Более того, при ρ = 1 имеем NS(J) ⊗ ℚ = ℚ·Θ и **Ker(NS(J) → NS(X)) = 0** — то есть у метода вообще нет исходного
материала: нужного нетривиального класса Z не существует, никакой квадратичной функции построить нельзя.
При ρ = 1 условие BD вырождается в классическое условие Шаботи–Коулмана r < g.

Проверка по другим полям (все безуспешны):
| поле F | ρ(J/F) | g + ρ − 1 | rank J(F) | вывод |
|---|---:|---:|---:|---|
| ℚ | 1 | 2 | 2 | не проходит |
| k = ℚ(√165) | 2 | 3 | 4 (= 2·rank E₁(k)) | не проходит |
| k(i) | 3 | 4 | 8 | не проходит |

Бюджет BBBM (1.4) над k: r + rk(𝒪_k^×) = 4 + 1 = 5 > 4 = [k:ℚ]·g — тоже не проходит.

**Это надо записать честно в FAMILY_SECTIONS.md: квадратичный Шаботи закрыт для всех шести сечений с w = +1.**

### 2bis.5. Единственная лазейка и как её проверить (дёшево)

Так как j(E₁) ∈ ℚ и j ≠ 0, 1728, кривая E₁ — **квадратичная закрутка** некоторой E₀/ℚ: E₁ = E₀ ⊗ δ, δ ∈ k*.
Тогда E₁^σ = E₀ ⊗ δ^σ = E₁ ⊗ (δ^σ/δ), и вся дихотомия сводится к одному square-class:

  **End⁰_ℚ(Res_{k/ℚ}E₁) ≠ ℚ  ⟺  E₁^σ ≅_k E₁  ⟺  δ^σ/δ ∈ k*²  ⟺  δ ∈ ℚ*·k*²  ⟺  E₁ опускается на ℚ.**

- Если **δ ∈ ℚ*·k*²** (E₁ — база-замена E₀′ = E₀ ⊗ δ с ℚ): тогда Res_{k/ℚ}(E₀′_k) ~ **E₀′ × E₀′^{(D)} над ℚ** — якобиан
  расщепляется, ρ(J/ℚ) = 2, кривая сечения **бисэллиптична над ℚ**, rank = r(E₀′) + r(E₀′^{(D)}) = 2.
  Тогда: (i) условие BD выполнено (2 < 3); (ii) **и вообще не нужен QC** — при ранге 1 + 1 работает готовый
  **Sage**-код Bianchi / Bianchi–Padurariu (`QC_bielliptic`), а при ранге 2 + 0 достаточно тривиального аргумента.
- Иначе (наш случай, δ = √165): ρ = 1, всё закрыто.

**Проверка на одну строчку кода для каждого сечения:** вычислить δ (отношение коэффициентов модели E₁ к модели E₀ над ℚ)
и проверить, лежит ли δ в ℚ*·k*². Эквивалентно: проверить, есть ли у кривой сечения **два ℚ-рациональных эллиптических
фактора** (бисэллиптичность над ℚ). По описанию проекта («бисэллиптична не над ℚ, а над k») ответ, по-видимому,
отрицательный для всех 11, но проверить стоит — это единственный сценарий, в котором сечение закрывается дёшево и в Sage.

Замечание: закручивать не поможет. Для d ∈ ℚ*: Hom_k((E₁⊗d)^σ, E₁⊗d) ≅ Hom_k(E₁^σ, E₁) = 0, то есть квадратичные
твисты кривой сечения над ℚ не меняют ни End⁰, ни ρ. «Правильный» твист (на √165) выводит E₁ в ℚ, но соответствующая
поверхность E₀^{(165)} × E₀ — уже не якобиан нашей кривой.

### 2bis.6. Вопрос 3: есть ли вариант QC не через NS?

**Короткий ответ: нет ничего запускаемого.**

- **BD II, Theorem 1.1** (обход через J ~ A^d × B, «Манин–Демьяненко»): условие
  ρ_f(A)·d + d(d−1)e(A)/2 − 1 > min{d(r − dim A), r² − dim(A)²}. У нас J ℚ-проста ⇒ d = 1, A = J, dim A = 2, r = 2, ρ_f = 1:
  слева 1 + 0 − 1 = 0, справа min{0, 0} = 0; 0 > 0 **ложно**. Не проходит.
- **Мнимое квадратичное умножение** (гипотетический вариант из вопроса координатора) — как показано в 2bis.3, даёт ρ = 1,
  то есть тоже не даёт QC. Отдельной ветки «QC для абелевых поверхностей с мнимым квадратичным умножением» **не существует**,
  и по структурной причине: метод BD работает ровно с Розати-симметричной частью End.
- **QC над числовым полем** (BBBM, arXiv:1910.04653) — см. таблицу в 2bis.4: ни над k, ни над k(i) бюджет не сходится.
- **Chabauty–Kim большей глубины** (Kim, Dogra, Hast, Corwin–Dan-Cohen) — теоретически покрывает случай ρ = 1,
  но **реализаций нет** (единственные вычислительно доведённые уровни — глубина 2, то есть QC).
- **Остаётся ровно то, что мы и делаем:** подъём на неразветвлённые накрытия (covering collections / two-cover descent),
  где ранг падает, плюс elliptic curve Chabauty на накрытиях и решето Морделла–Вейля. См. п. 2bis.7.

### 2bis.7. Вопрос 4: схема Bars–González–Xarles при rank = [k:ℚ]

**F. Bars, J. González, X. Xarles, «Hyperelliptic parametrizations of Q-curves», Ramanujan J. 56 (2021), 103–120;
препринт arXiv:1910.10545.** Код (Magma): github.com/XavierXarles/HyperellipticParametrizationsQcurves (`SuccesfulCurves.m`).
Кривые: X₀*(N) = X₀(N)/B(N) (полная группа Аткина–Лемера), N бесквадратно; при бесквадратном N гиперэллиптичность ⟺ род 2.
Мотивация — теорема Элкиса: всякая ℚ-кривая без CM изогенна над ℚ̄ ℚ-кривой, отвечающей рациональной точке некоторой X₀*(N).

**Метод — дословно (§3.2):**
> «In order to determine the rational points of the curves X₀*(N) we will use the so-called **elliptic Chabauty method,
> which uses a Chabauty procedure on a finite set of unramified 2-coverings of the curve**.»
> «first one computes the finite set of twists C_ξ of the unramified coverings of the curve X with Galois group ≅ (ℤ/2ℤ)^{2g}
> which have points locally for any prime p; this is completely analogous to the 2 descent for elliptic curves, as described in [BrSt09]…
> Now, the jacobian of any of this curves has quotients isomorphic to the **Weil restriction of elliptic curves E_ξ defined over some
> number fields K**, and the rational points in C_ξ(ℚ) give points in E_ξ(K) whose image with respect to a given map φ_ξ: E_ξ → **P**¹
> is in **P**¹(ℚ); this is the necessary data for the elliptic Chabauty function, which computes the set of points in E_ξ(K) verifying
> this condition **if rank_ℤ(E_ξ(K)) < deg(K/ℚ). In practice, the fields K we need are the minimal field of definition of some fixed
> factorization f(x) = g(x)h(x) where g(x) has degree 4.**»

Решето Морделла–Вейля — только запасной вариант и **не понадобилось**: «For our curves this never happens».

#### ⚠️ Главное для нас: их rank и наш rank — это РАЗНЫЕ величины

У BGX E_ξ — эллиптический фактор якобиана **неразветвлённого 2-накрытия** C_ξ, а K — поле определения разложения
f = g·h с deg g = 4 (его можно выбирать!). Наша пара (E₁, k) — это уровень **самой кривой сечения**, до спуска.
**То, что rank E₁(k) = 2 = [k:ℚ], ничего не говорит о том, пройдёт ли схема BGX на 2-накрытиях.**
Иначе говоря, наш вывод 2026-09-11 «эллиптический Чабо неприменим при ранге 2» верен только для прямого захода
на E₁ над k, а не для схемы «covering collection + elliptic Chabauty».

Наглядно: для C₁₆₁ (= X₀*(161), та самая кривая из BDMTV §5.4) rank Jac(ℚ) = 2 = g, классический Шаботи мёртв,
а BGX берут **кубическое** K = ℚ[x]/(x³ − 2x² + 3x − 1), и там rank E_ξ(K) = 2 < 3 = [K:ℚ] — условие спасено
**ростом степени поля**, а не падением ранга. (Из их лога: `Rank for twist= 1 = [ 2, 2 ]` — ранг ровно 2.)

#### Три рабочих рецепта при rank = степень

**A. Сменить K (самый дешёвый).** Перебрать все разложения f = g·h с deg g = 4 над разными подполями поля разложения f.
BGX, X₀*(85) — буквально наш сценарий:
> «The twists corresponding to the points with x-coordinate 3/2 and −4/3 have rank 1, and Chabauty method succeeds.
> **But the ones corresponding to the points with x-coordinate 0 and 2 have rank 2.** If we adjoint a root of x⁴−2x³+3x²−6x+5
> we get a field K₄ where f(x) has 4 roots and a degree 2 factor… **we finally get that the remaining twist have jacobian of rank 1…
> and Chabauty computations succeeds.**»
То есть rank 2 = [K₂:ℚ] → провал; переход к K₄ степени 4 → rank 1 < 4 → успех.
Тот же приём у Брюина (Crelle 562, §3): «for C₇ and C₉, there is another choice that would yield elliptic curves over a cubic
extension. Those elliptic curves have rank 3, which is too high for the method to work.»

**B. Подняться на неразветвлённое 2-накрытие и искать K заново.** BGX, X₀*(390):
> «we tried fields of degree 1 and 2 and the Chabauty condition was not fulfilled, and we had to go to a degree 8 extension,
> where the rank computation did not succeed. Instead… **The curve X′ determined by these equations is an hyperelliptic curve of
> genus 3**… **Now we apply the above method to this new curve X′**… **Over K = ℚ(√5)** the defining hyperelliptic polynomial factors
> as a product of two degree 4 polynomials. **For every twist one of the two quotient elliptic curves of the corresponding covering
> has rank one, and the Chabauty computation succeeds.**»
Заметим: сработало именно **вещественное квадратичное** поле, на уровне накрытия.
Ещё ближе к нам — **Flynn–Wetherell, «Covering collections and a challenge problem of Serre», Acta Arith. 98 (2001), 197–205**
(X⁴+Y⁴=17Z⁴): род-2 кривая с «the Jacobian has rank 2 over ℚ, and so the problem **just barely eludes this attack**» (rank = g = 2!)
закрывается covering collection → кривая рода 1 над K = ℚ(√2, √34), [K:ℚ] = 4, rank 1 < 4.

**Важное предостережение** — N. Bruin, E. V. Flynn, «Towers of 2-covers of hyperelliptic curves», Trans. AMS 357 (2005), 4329–4347:
> «Suppose we have a curve C over a number field K with rk Jac(C)(K) ≥ genus(C) and a cover D such that Jac(D) is isogenous over K
> to Jac(C) ⊕ E₁ ⊕ … ⊕ E_n… **If all the E_i have positive Mordell–Weil rank, then the cover D/C is essentially useless from the point
> of view of Chabauty's method.** One encounters this situation if one considers a hyperelliptic curve C of genus 2 with Jac(C)[2]
> pointwise defined over K and D a pullback of an embedding of C in Jac(C) along multiplication by two.»
**Это ровно наша конфигурация** (род 2, полное 2-кручение над k) — поэтому одношаговое 2-накрытие может не помочь,
и Брюин–Флинн предлагают **башню** (повторный спуск через другое подкрытие C′).

Важное уточнение общей теории (Bruin, диссертация Leiden 1999 / CWI Tract 133, §5.1), дословно:
> «**If rk(J(K)) is high, then a different method is needed. If we have a finite number of covers D → C over K such that the images
> of D(K) cover C(K), then we can try to determine D(K) through the same process. In general, genus(D) > genus(C) and
> rk(Jac(D)(K)) ≥ rk(Jac(C)(K)). Therefore, it is possible Chabauty methods are applicable to the D.**»
**Ранг на накрытии не падает — он ≥ исходного.** Выигрыш чисто относительный: растёт род, якобиан дробится,
эллиптические факторы живут над бо́льшими полями, и условие rank < [K:ℚ] может выполниться.

**C. Quadratic Chabauty вместо elliptic Chabauty** — для нас закрыт (ρ = 1, см. 2bis.3–2bis.4).
Именно этим путём BDMTV решили C₁₆₁ (у них ρ(J) = 2, RM), нам он недоступен.
Adžaga–Chidambaram–Keller–Padurariu (arXiv:2203.05541, продолжение BGX на все 64 гиперэллиптических X₀*(N)) формулируют
развилку так:
> «**Note that the genus of the coverings will grow, so the Chabauty condition is more likely to be satisfied at the expense of
> increasing the degree.**… **One usually suspects this method to work if the degree of K is large because the ranks of elliptic curves
> are expected to be small**… **However, elliptic curve Chabauty often fails for algorithmic reasons.** Namely, if the étale algebra K
> has too large a degree, to perform the 2-descent to compute A(**Q**), one would need to compute the class group of a number field of
> large degree. Such a computation is infeasible, often even assuming GRH. **In these cases we can use the quadratic Chabauty method.**»
(Они же отмечают опечатку у BGX: бесквадратных уровней не 36, а 39 — пропущены 166, 255, 330.)

#### Точная формулировка условия Брюина

**N. Bruin, «Chabauty methods using elliptic curves», Crelle 562 (2003), 27–49**
(препринт: cecm.sfu.ca/~nbruin/bruin_crelle.pdf). Условие **не оформлено как теорема**; §4.2, дословно:
> «Let ℚ ⊆ K ⊆ L be number fields and let φ: E → **P**¹ be an elliptic cover defined over L. In this section we propose a method for
> determining the L-rational points G on E such that φ(G) is K-rational… **The method we explain here might give a sharp bound on the
> number of such G if rk(E(L)) < [L : K].**»
Обратите внимание: условие **относительное** — rk E(L) < [L : K] (у BGX K = ℚ, отсюда их форма).
Механизм (Lemma 4.1): нужна матрица m × r ранга r, то есть **r < m ≤ [L:K]**; при r = [L:K] метод структурно ломается.
В диссертации то же: «**A necessary condition for the method to work is that rk(E(K)) < [K:ℚ].**»

Отдельно: **Flynn–Wetherell, «Finding rational points on bielliptic genus 2 curves», manuscripta math. 100 (1999), 519–533** —
условие (2.12) «Suppose that the rank r of E(ℚ(α)) is less than d, the degree of ℚ(α):ℚ»; и честный пример провала:
> «**We found precisely one curve where the method failed**, namely: C: Y² = (X²+1)(X²+3)(X²+7)… **Hence we cannot reduce our
> eligible x = X² down to a finite number of possibilities, and the method fails.**»
Расширение: **Mourao, «Extending elliptic curve Chabauty to higher genus curves», manuscripta math. (2013), arXiv:1111.5506** —
ослабляет условие до rank ≤ dg − 1; при g = 1 это то же самое rank < d, **выигрыша на эллиптическом уровне нет**,
выигрыш только на накрытиях большего рода.

#### Софт

Всё ядро — **Magma**: `TwoCoverDescent` (Bruin–Stoll), `Chabauty(mwmap, u : IndexBound:=...)`, `Chabauty0`,
`DescentInformation(Em: RankOnly:=true)`, `Saturation`, `IsPSaturated`, `EvaluateByPowerSeries`.
Ключевая строка кода BGX — буквальная проверка условия:
```magma
rank, gs := DescentInformation(Em: RankOnly:=true);
if rank[1] lt Degree(K) then     // rank E(K) < [K:Q]
```
Если условие нарушено, твист **молча пропускается** — надо пробовать другое K или другую кривую.
В Sage `Chabauty`/`TwoCoverDescent` нет; D. Hast (arXiv:2009.10313, github.com/HastD/twocover-descent) прямо пишет,
что чисто-Sage реализация неполна, ядро вызывает Magma. В Sage/PARI доступно только `E.rank_bounds()` /
`simon_two_descent()` над числовыми полями — **этого хватает на разведку**: перебрать кандидатов K и посмотреть,
где ранг падает ниже степени.

#### Практический план для наших шести сечений (w = +1)

1. Взять модель кривой сечения C: y² = f(x) над ℚ и запустить `TwoCoverDescent` — получить твисты δ.
2. Для каждого δ перебрать **все** поля K = поле определения разложения f = g·h с deg g = 4, включая **кубические
   и квартичные** подполя поля разложения f — а не только наше k = ℚ(√D). Проверять буквально `rank < Degree(K)`.
   **Разведка возможна в Sage** (`rank_bounds` над K), решающий шаг — Magma.
3. Если ни одно K не проходит — подняться на накрытие рода 3 (рецепт B, как X₀*(390)) и повторить;
   держать в уме предупреждение Брюина–Флинна про «essentially useless» накрытия при полном 2-кручении и, если надо, строить башню.
4. Quadratic Chabauty как план C — **у нас закрыт** (ρ = 1).

Готовый шаблон: функция `MyEllipticChabauty(g, deltas, uXrat)` в `SuccesfulCurves.m`
(github.com/XavierXarles/HyperellipticParametrizationsQcurves), там же логи всех 19 случаев, включая N = 161.
Это задача для Codex (Magma).

---

## 3. Поверхность магического квадрата: Auel–Singer

**A. Auel, B. Singer, «The algebraic geometry of 3-by-3 magic squares of squares», arXiv:2609.09351 (8 сентября 2026), 43 стр.**
MSC 14G05, 14J29, 14J28, 14Q10, 14J50. Полный текст доступен в HTML-версии arXiv (https://arxiv.org/html/2609.09351v1).

Обозначения статьи: V ⊂ P⁸ — поверхность (полное пересечение шести квадрик, степень 64, 256 узлов A₁), Ṽ — разрешение,
U ⊂ V — «distinct locus» (все девять записей различны).

### 3.1. Нерегулярность и инварианты — **q = 0**

§2, цитата:
> «For any complete intersection X of dimension d and for any 0 < i < d, we have that h^i(X, 𝒪_X) = 0, so
> q(Ṽ^an) = h^{0,1}(Ṽ^an) = h^{1,0}(Ṽ^an) = 0.»

**Theorem 2 (Hodge Numbers, Betti Numbers, and Chern Numbers)**, цитата:
> «The complex manifold Ṽ^an is simply connected and has Hodge diamond 1 / 0 0 / 111 544 111 / 0 0 / 1.
> This gives Betti numbers b₀ = b₄ = 1, b₁ = b₃ = 0, b₂ = 766. Furthermore, the singular cohomology of Ṽ^an is torsion-free,
> and the Chern numbers are c₁²(Ṽ^an) = 576 and c₂(Ṽ^an) = 768. The lattice H²(Ṽ^an, ℤ) with the intersection pairing is odd and
> unimodular of signature (223,543), hence isomorphic to ⟨1⟩^⊕223 ⊕ ⟨−1⟩^⊕543.»

Итого: **q = 0**, p_g = 111, h^{1,1} = 544, K² = 576, c₂ = 768; χ(𝒪) = 112 (в статье явно не выписано, пересчёт по Нётеру).
**Lemma 1.1**: «The resolution Ṽ is of general type» (ω_Ṽ ≅ π*𝒪_V(3), разрешение узлов A₁ крепантно).

Важное следствие (конец §2), цитата:
> «Finally, we remark that since q(Ṽ) = 0, the Albanese variety of Ṽ, Alb(Ṽ), is a point. Hence, Ṽ admits no interesting morphisms
> to abelian varieties.»

**Theorem 3 (Picard Rank)**, цитата: «The resolution Ṽ has torsion-free geometric Picard group, with geometric Picard rank
518 ≤ ρ(Ṽ_Q̄) ≤ 544.» (Нижняя граница — явная решётка из 1204 дивизоров; Remark 5.2: авторы ожидают максимальный ранг 544,
но «The methods of Sections 3 and 4 appear to be saturated… producing the remaining classes [will] require Hodge-theoretic methods».)

### 3.2. Кривые малого рода; кривые степени 8

**Prop. A.1 (Appendix A), цитата:** «The Fano scheme of lines F(V) of V, the variety of 3×3 magic squares of squares, is empty.»
**Прямых на V нет.**

**Proposition 3.1 (вырожденный локус), цитата:**
> «The nondistinct locus Z ⊂ V decomposes into Z₃ ∪ Z₅ ∪ Z₇, where each component Z_i generically corresponds to squares with i distinct entries.
> All irreducible components of Z are smooth geometrically connected curves defined over ℚ…»
с таблицей: **Z₃ — 128 компонент, род 0, степень 2; Z₅ — 32 компоненты, род 5, степень 8; Z₇ — 16 компонент, род 49, степень 32.**
Итого 176 неприводимых компонент. 128 коник Z₃ — это «квадраты в арифметической прогрессии» B² + C² = 2A².

**Proposition 3.3 (эллиптические кривые), цитата:** «The M = 0 hyperplane section decomposes, over ℚ(i), into the union of
16 smooth elliptic curves all isomorphic to y² = x³ − x.» (Prop. 3.6: Z_M — 16 компонент, **род 1, степень 4**.)

**Кривые степени 8 — все рода 5.** В статье они встречаются несколько раз: Z₅ (32 компоненты, род 5, степень 8, Prop 3.1);
Z_middle (32 компоненты, род 5, степень 8, Prop 3.6); Z_{B,C} над ℚ(i,√3) и Z_{B,M} над ℚ(√−2) (по 64 компоненты, род 5, степень 8, Prop 3.7);
компоненты 384 дивизорных классов из §4.1 — «smooth curves of genus 17 and degree 16, or genus 5 and degree 8».
**Ни одной кривой степени 8 рода 0 или 1 в статье нет**; кривые рода 1 имеют степень 4, рода 0 — степень 2.
(Замечание: «degree 8» в §4.2 — про другое: «magic octic K3 surfaces», K3 степени 8 в P⁵.)

Классификация кривых малого рода **не завершена** — Remark A.2, цитата:
> «This trivially implies that no lines meet the distinct locus U. We record it because the question of low-degree rational curves on U —
> rather than on all of V — is what is directly relevant to the Diophantine problem. By contrast, V contains many conics.
> The 128 rational conics arising from the 3-distinct-entry locus Z₃ in Section 3.1 are all contained in the closed locus Z₃ ⊂ V∖U.
> We leave open the question of whether these are all the conics contained in V; we know there are finitely many by [12].»

### 3.3. Ланг / Бомбьери–Ланг; конечность исключительного множества

Собственного результата о конечности исключительного множества Ланга **в статье нет**. Есть два места.

Введение, цитата:
> «The resolution Ṽ is of general type, so the sparsity of 3×3 magic squares of squares is consistent with the Bombieri–Lang conjecture…
> The 3×3 case is the only case where this happens; the variety associated to 4×4 magic squares of squares is a Calabi–Yau sevenfold and the
> varieties associated to n×n squares of squares for n ≥ 5 are Fano, so they should contain many rational points. … Furthermore, Bruin, Thomas,
> and Várilly-Alvarado [12] have shown Ṽ to be algebraically quasi-hyperbolic, i.e., it contains finitely many curves of genus 0 and 1.
> Since these are the only curves that could have infinitely many rational points by Faltings' theorem, it confirms that rational points on V
> are relatively sparse.»

([12] = N. Bruin, B. Thomas, A. Várilly-Alvarado, «Explicit computation of symmetric differentials and its application to quasihyperbolicity»,
Algebra Number Theory 16 (2022), 1377–1405 — **это и есть доказанный результат о конечности числа кривых рода 0 и 1**, не Auel–Singer.)

После Lemma 1.1, цитата (важное предупреждение):
> «The Bombieri–Lang conjecture [30] predicts that the rational points of Ṽ and V are contained in a Zariski closed subset, i.e., a finite union
> of curves and points. This would guarantee some open subset V′ ⊂ V with no rational points. However, this is still not enough to solve the
> 3×3 magic square of squares problem; since a square having distinct and nonzero entries is an open condition, solving the problem is equivalent
> to showing that a specific open subset U ⊂ V has no rational points.»

### 3.4. Семейства с бесконечным числом точек

**Единственный доказанный случай бесконечности — Proposition 3.4, цитата:**
> «For d a congruent number, the distinct locus U has infinitely many points with M = 0 over ℚ(i, √d).»

То есть бесконечные серии есть только **при M = 0** (нулевой центр) и только над биквадратичными полями ℚ(i,√d) — над ℚ бесполезны.
**Remark 3.5**: «the elliptic curve y² = x³ − x has Mordell–Weil group structure ℤ/2ℤ ⊕ ℤ/4ℤ over ℚ(i)» ⇒ над ℤ[i] решений с M = 0 нет.

**Эллиптических расслоений на самой Ṽ в статье нет** (единственное упоминание elliptic fibration — у Бремнера, на K3-поверхности орбиты 7).
Рангов Морделла–Вейля пучков на V авторы не считают.

Известные «почти-магические» семейства, как они изложены в статье:
- **Lucas (1891)** — параметризация всех 3×3 магических квадратов (формула (1) во введении); он же доказал невозможность бимагического 3×3.
- **Bremner [10]** — квадрат квадратов над ℚ(√3,√133) с суммой 1596 («a magic square of squares over the totally real field of smallest degree known»);
  семейство с нулевым центром над ℚ(i,√(q³−q)).
- **Новые семейства авторов** (§3.2, M = 0): через пифагоровы тройки — над ℚ(i,√(2(m⁴+n⁴))) (пример ℚ(i,√34)) и над ℚ(i,√(m⁴+4n⁴))
  (пример ℚ(i,√5), «this discriminant is smaller than is achievable for Bremner's family»).
- **5 квадратов** — §4.1: проекции V → P⁴, 23 орбиты, **5 геометрических классов магических квартических дель Пеццо** (Prop 4.2, Table 1).
  Дель Пеццо рациональны ⇒ квадратов с 5 квадратами бесконечно много; Remark 4.5 приводит подъём точек высоты 28.
- **6 квадратов** — §4.2: 16 орбит, «magic octic K3 surfaces» (полные пересечения трёх квадрик в P⁵); Prop 4.11: «Up to automorphism,
  there are 14 magic octic K3 surfaces»; Prop 4.8: ρ(X̃_Q̄) = 19 для орбиты 7; у гладкой поверхности Бремнера (орбита 3) ρ = 20.
- **7 квадратов** — известный единственный пример (центр 425²) лежит в слое над рациональной точкой [373:289:565:425:205:527] K3-поверхности X;
  цитата: «and it's purely a coincidence that E is a square in this fiber. Previous work of Bremner [11] also shows that squares of 7 squares
  arising from this orbit are rational points on high genus hyperelliptic curves on X, giving an alternate explanation.»
  ⇒ **бесконечного семейства с 7 квадратами геометрия не обещает**.

### 3.5. Что геометрия отсекает — и что это значит для нас

- **Вырожденный локус описан исчерпывающе** (Prop 3.1, Prop 3.6, Prop 3.7): 176 + 64 компонент, и только Z₃ (128 коник, род 0)
  и Z_M (16 кривых рода 1) могут нести бесконечно много точек — **и обе серии вырождены** (3 различные записи; M = 0).
  Remark 3.2, цитата: «The fact that the components of Z₅ and Z₇ are of high genus is consistent, via Faltings theorem, with older results
  to the effect that there are no integer 3×3 magic squares of squares with 5 or 7 distinct entries.»
- **На U (нужный нам открытый локус) не найдено ни одной кривой рода 0 или 1**; по [12] их заведомо конечное число.
  Это значит: **наш поиск по одномерным семействам (сечения) на U не может быть закрыт «одной кривой рода 0 с бесконечным числом точек»** —
  таких кривых там, по-видимому, нет; наши сечения — это кривые рода ≥ 2 (у нас род 2), и Фальтингс/Шаботи применимы по существу.
  Это **подтверждает выбранную стратегию** (конечность точек на каждом сечении + перебор сечений), а не отсекает её.
- **Ограничение методов**: Ṽ односвязна (Theorem 2) ⇒ цитата из введения:
  «Since Ṽ is simply connected by Theorem 2, it admits no nontrivial connected finite étale covers. Thus any étale input for a potential
  étale Brauer–Manin obstruction on U must come from covers ramified along its complement.»
  Плюс Alb(Ṽ) = точка ⇒ никаких нетривиальных отображений в абелевы многообразия. **Глобальных обструкций «сверху» ждать не приходится**;
  работать надо с разветвлёнными накрытиями — что мы и делаем (2-накрытия, Richelot).
- §3.3 фиксирует тупик в одном направлении: для A = ±√d·M при d > 3 «these hyperplane sections fail to split further and become
  linearly dependent; their intersection matrix is rank 1.»

### 3.6. Чего в статье Auel–Singer нет

1. χ(𝒪_Ṽ) явно не выписана; плюригенусы, каноническая модель не обсуждаются.
2. Нет описания/оценки исключительного множества Ланга — только ссылка на гипотезу Бомбьери–Ланга и на [12].
3. Нет полной классификации кривых рода 0 и 1 на V (открытый вопрос: все ли коники — это 128 коник Z₃).
4. Нет эллиптических расслоений на Ṽ и расчётов рангов Морделла–Вейля.
5. Нет утверждений «кривых степени 8 рода ≤ 1 не существует» (доказано только отсутствие прямых).
6. Нет результатов над ℚ, дающих бесконечные семейства.

---

## Чего найти не удалось (сводно)

**По вопросу 1:**
- **Полный текст Cassels 1998** (Crelle 494, 101–127) — нигде в открытом доступе; DigiZeitschriften закрыт с 31.12.2025,
  GDZ том не отдаёт, препринта нет. Точные номера лемм известны только по цитированиям: Lemma 7.2 (локальное зануление),
  Lemma 7.4 (определение спаривания), §2 (модели 2-накрытий). Формулировка выше взята из двух независимых пересказов,
  которые между собой согласуются (Shukla–Stoll над числовым полем; Zhang над ℚ с полным 2-кручением) плюс третья,
  структурно эквивалентная форма (Fisher–Newton / Visse).
- Не проверено, работает ли Magma `CasselsTatePairing` **над числовым полем** (нам нужно k = ℚ(√165)).

**По вопросу 2:**
- Обзор Balakrishnan «Quadratic Chabauty» в трудах **ICM 2022** — существование не подтверждено; вероятно имелся в виду конспект
  Balakrishnan–Müller, «Computational tools for quadratic Chabauty» (Arizona Winter School 2020),
  math.bu.edu/people/jbala/2020BalakrishnanMuellerNotes.pdf (сорвался fetch по SSL — точных цитат нет).
- **arXiv:1906.04143** с содержанием «Quadratic Chabauty for (bi)elliptic curves» — такой работы нет.
  Нужные: Bianchi arXiv:1904.04622 и BBBM arXiv:1910.04653.
- **Ни одной работы, где квадратичный Шаботи применяется к кривой рода 2 с J ~ Res_{k/ℚ}E при вещественном квадратичном k.**
- Публично выложенного репозитория к BBBM (QC над числовым полем, Example 7.2, ℚ(√34)) найти не удалось.
- Verbatim формулировки ограничения `rank ≥ degree` из Magma handbook — страница отдала 401; цитата взята из сниппета поиска.
- Точной теоремы «End⁰(Res_{k/ℚ}E) = ℚ ⟺ E не ℚ-кривая» с номером — не найдено; опора косвенная (Quer, Ribet).

**По уточнению 2bis:**
- **Полного текста Quer, PLMS 81 (2000), 285–317** получить не удалось (Oxford Academic за подпиской); вывод про
  End⁰(Res) опирается на Ribet §6–7 (свободный текст) и González–Guàrdia–Rotger arXiv:math/0409352 (Lemma 3.6, Def. 2.7).
- Формула ρ(A/ℚ) = dim(Розати-симметричная часть End⁰_ℚ(A)) — классика (Mumford §20–21; Birkenhake–Lange Prop. 5.2.1),
  но дословной цитаты «над ℚ, Галуа-эквивариантно» я не выписывал — это стандартное следствие ℚ-рациональности Θ.
- **Журнальной версии BGX (Ramanujan J.)** — только за paywall; нумерация предложений в журнале может отличаться от arXiv:1910.10545
  (там Proposition 1 в §2, Proposition 2 в §3.2). **Перед цитированием в публикации сверить нумерацию.**
- Страница Magma handbook «Elliptic Curve Chabauty» отдаёт 401/403 — сигнатуры `Chabauty(MWmap, Ecov)` взяты из кода BGX.
- Отдельного PDF CWI Tract 133 (2002) нет; доступна только диссертация Leiden 1999 (cecm.sfu.ca/~nbruin/thesis.pdf).
- **Ни одной работы, где написано буквально «rank = degree, поэтому мы применили решето Морделла–Вейля»** — решето в этих
  схемах используется для доказательства *пустоты* C_ξ(ℚ), а не для обхода ранга.
- Все ранги и логи BGX приведены по опубликованным логам авторов (Magma в окружении нет) — статус **numerical**, не перепроверено.

**По вопросу 3:**
- Вторичных источников (обсуждения, цитирующие работы) на Auel–Singer нет — статья от 08.09.2026.
- В самой статье нет: χ(𝒪) явно, описания/оценки исключительного множества Ланга, полной классификации кривых рода 0 и 1,
  эллиптических расслоений на Ṽ, утверждений о кривых степени 8 малого рода, семейств с бесконечным числом точек над ℚ.
</content>
