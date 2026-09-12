# Поиск (Grok), 12.09.2026

Ответ на `QUESTION_TO_GROK_2026-09-12_SEARCH.md`, пункт 1 (оставлен открытым после Gemini) плюс короткие уточнения к 2–4.
Не расчёты. Sage не менял. Сайт не трогал.

Две классификации **не смешивать** (Gemini это уже поправил, здесь закреплено):
- Boyer **7.I–7.VIII** — восемь конфигураций **семи** квадратных клеток (Search.pdf, fig. 5).
- Bremner **I–XVI** — шестнадцать конфигураций **шести** квадратных клеток (Acta Arith. 99 (2001), Fig. 1).

---

## 1. Кто что перебирал (покрытие)

### 1.1 Chaffin — семиклеточные Boyer, не шестёрки Bremner

Источник: https://benchaffin.com/magic-squares/magic-squares.html (работа ~2008, страница 2019).

**Что именно.** Все конфигурации Boyer **7.I–7.VI** (квадратный центр) и **7.VII–7.VIII** (неквадратный центр). Метод: генерирует числа `x` с простыми только 4k+1, раскладывает `2x` и `2x²` в суммы двух квадратов, расставляет пары в клетки по схеме 7.x.

**Границы — на центральную клетку E, не на параметры семейства и не на все девять корней:**
> «There are no more magic squares which use 7 or more squares which have a **non-square central cell up to 10¹⁴**, or a **square central cell up to 10²⁸**.»
> «The search considered over 3 trillion values for E, and checked about 875 trillion magic squares.»

Для квадратного центра: `E = x²`, нужно ≥ 3 разложения `2x²` (одно из них `x²+x²` отбрасывается). Для неквадратного: `E = x`.

**Вывод для нас.** Это закрывает **достройки шестиклеточных до семи** с данным центром в этих границах. Это **не** перебор 16 однопараметрических семейств Bremner по высоте параметра λ. Наш MW-перебор по G-семействам и Chaffin — разные оси: он режет по величине центра, мы — по высоте на кривой фиксированного сечения.

### 1.2 Buell — только песочные часы, и только взаимно простые линии

Duncan A. Buell, *A search for a magic hourglass*, preprint 1999 (ссылка с Boyer: http://www.multimagie.com — PDF «A search for a magic hourglass»). Цитата через Pierrat–Thiriet–Zimmermann 2015:

> «Duncan Buell shows … if a solution exists, its center cell is larger than 25 · 10²⁴.»
> «he assumes that in each of the two diagonals and the central column, the three entries are coprime, which does not necessarily hold.»

Hourglass = три линии через центр = 7 квадратов (Boyer 7.x типа «три AP через центр»). **Не** все 16 шестёрок.

Pierrat–Thiriet–Zimmermann, *Magic squares of squares* (2015), https://members.loria.fr/PZimmermann/papers/squares.pdf — сняли coprimality, нашли hourglass с центром из 10 цифр, у Buell/Pech modulo пусто до 10¹²–10¹³.

**Вывод.** Граница 25·10²⁴ **не** покрывает общий случай и **не** покрывает Bremner I–XVI.

### 1.3 Morgenstern / Boyer — семиклеточный метод 2010, не каталог I–XVI

Morgenstern, окт. 2010, http://www.multimagie.com/English/Morgenstern08.htm: **полный** метод для семи квадратов через **шесть** конфигураций из пяти квадратов (две AP с общей клеткой). Конфиги 1–3 = квадратный центр = наши G1/G2/G3; 4–6 = неквадратный центр. Это полнота постановки, не таблица «конфиг I Bremner перебран до …».

Публичные стены (сводка uwe-schwarz independent-2026-audit, 2026-08-22, по первоисточникам): Morgenstern entry bound **10¹⁴** и множество квадратов mod 2⁹⁰. Boyer ~2,5·10²⁷ по центру (из NOVELTY_CHECK). Отдельных границ **по каждому из I–XVI** в открытом доступе **не нашёл**.

### 1.4 Bremner 2001 сам

Для **каждой** из I–XVI: бесконечные семейства шести квадратов (эллиптическое расслоение); седьмая клетка = гиперэллиптическая f(t)=□. Компьютерный поиск по этим кривым дал **только** известный семиклеточный пример. Это поиск рациональных точек на кривых, не brute по центру.

Конфиг **VI** = полный крест + один угол (уже в `GROK_INFO_7_8_SQUARES.md`).

### 1.5 После 2019 (кроме Auel–Singer)

| кто | что | граница | конфиги |
|---|---|---|---|
| Wesolowski, dec 2019 / jan 2020, Boyer Search page | одно 6-квадратное семейство (сужение G1) | odd n ≤ 2·10⁵; even n — **доказано**, красные клетки не квадраты | не 16 конфигов, одно семейство |
| Wesolowski, май 2020 | второе 6-квадратное семейство | красные клетки не квадраты для n < 10¹⁰ (сумма < 9.72·10⁸²) | то же |
| Chaffin | см. 1.1 | 10²⁸ / 10¹⁴ на E | 7.I–7.VIII |
| uwe-schwarz / mystimath 2026 | 9/9 корень ≤ 10⁶; p²qr классы | не 7-клеточные конфиги Bremner | центр |
| philthompson.me 2025–02.2026 | полный 9/9 | корни ≤ 512 | все 9 клеток |
| Pierrat–Thiriet–Zimmermann 2015 | hourglass без coprime | примеры, не запрет | 7 клеток, часы |
| Auel–Singer 2026 | геометрия поверхности | нет перебора 16 конфигов | — |

**Итог пункта 1.** Публичного «конфиг Bremner k закрыт до высоты H» **нет**. Есть: (а) семиклеточные Boyer 7.x по величине центра (Chaffin); (б) hourglass с оговоркой coprime (Buell); (в) кривые Bremner 2001; (г) два семейства Wesolowski. Наш brute 13,8 млрд по 16 конфигам **не дублирует** Chaffin один в один, но для **семиклеточной** цели с квадратным центром Chaffin уже жёстче по E. Имеет смысл не расширять box-scan I–XVI, а арифметику G-сечений (уже начата).

---

## 2–4. Уточнения к ответу Gemini (коротко)

**Magma calc.** Согласен: ~60 с, `TwoCoverDescent` / род-2 `Chabauty` там не живут. Открытой реализации two-cover descent вне Magma по-прежнему не нашёл (Bruin–Stoll в Magma; mdmagma — эллиптический Чабо). Письмо Stoll/Bruin — правильный канал, не онлайн-калькулятор.

**CTP.** PARI `casselspairing` только над ℚ (Allomber). Magma `CasselsTatePairing` над полями — Donnelly, для 2-накрытий эллиптических. Наша открытая реализация над ℚ(√D) с полным 2-кручением — ниша (уже в NOVELTY_CHECK).

**Auel–Singer.** arXiv:2609.09351, 8.09.2026, препринт. «Comments welcome!». Ṽ односвязна ⇒ этальный Brauer–Manin недоступен (не повторять ошибку Gemini). Контакты — `CONTACTS_2026-09-12.md`. Boyer: mailto живой на сайте, апдейт страницы с июня 2020.

**Wesolowski even n.** Boyer, SquaresOfSquaresSearch: «When n is an even positive number, we can prove that the cells in red can't be squared integers.» Это **доказательство** на его семействе, не поиск. Текст доказательства на странице не выписан — только факт. Имеет смысл восстановить (2-адика / mod 4) и перенести на объёмлющее G1.

Конец.
