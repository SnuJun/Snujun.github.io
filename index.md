# Welcome

Hello, I’m **Junwoo Lee**.  
I recently received my M.S. in Statistics from Seoul National University and will begin my Ph.D. in Statistics at Virginia Tech in September 2026.  
I am currently working as a research assistant in the Multivariate Statistics Lab.

---

### Education

- **Virginia Tech**  
  *Ph.D. in Statistics (Incoming)*  
  *Sep 2026 – Present*

- **Seoul National University**, College of Natural Sciences  
  *M.S. in Statistics*  
  *Mar 2024 – Feb 2026*  
  Advisor: Johan Lim  
  Focus: High-dimensional statistical inference and covariance estimation  

- **Konkuk University**, College of Social Sciences  
  *B.A. in Applied Statistics*  
  *Mar 2017 – Feb 2024*

- **Republic of Korea Army**  
  *Military Service*  
  *Mar 2019 – Mar 2021*

---

### Projects & Activities

#### Undergraduate
- [Machine Learning Club *Kuggle*](/kuggle/){:target="_blank" rel="noopener"}
- [Economc Data Analysis](/assets/pdf/economic_data_final_presentation.pdf){:target="_blank" rel="noopener"}
- [Regression Project](/assets/pdf/regression_project_final_report.pdf){:target="_blank" rel="noopener"}






#### Graduate
- [Linear Shrinkage for Positive Definite Covariance Matrix Estimation](https://github.com/SnuJun/Snujun.github.io/tree/main/assets/LPD)
- [Deep Learning Project: Methods for Interpreting Deep Learning Models](/assets/pdf/Methods%20for%20Interpreting%20DL%20models.pdf){:target="_blank" rel="noopener"}



---

### Links
- [About](/about)  
- [RSS Feed](/feed.xml)

## Posts
{% for post in paginator.posts %}
- **[{{ post.title }}]({{ post.url }})** <small>({{ post.date | date: "%Y-%m-%d" }})</small>
  {% if post.excerpt %}<br>{{ post.excerpt | strip_html | truncate: 120 }}{% endif %}
{% endfor %}

{% if paginator.total_pages > 1 %}
<nav>
  {% if paginator.previous_page %}[← Newer]({{ paginator.previous_page_path }}){% endif %}
  <span> Page {{ paginator.page }} / {{ paginator.total_pages }} </span>
  {% if paginator.next_page %}[Older →]({{ paginator.next_page_path }}){% endif %}
</nav>
{% endif %}
