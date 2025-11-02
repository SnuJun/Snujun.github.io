# Welcome

Hello, I’m **Junwoo Lee**.  
I am currently pursuing a Master’s degree in Statistics at Seoul National University,  
where I study **high-dimensional statistics** and related topics in modern inference.

---

### Projects & Activities

#### Undergraduate
- [Regression Project](https://github.com/SnuJun/Regression-Project)  
- [Economic Data Analysis](https://github.com/SnuJun/Economic-Data-Analysis)  
- [Machine Learning Club *Kuggle*](https://github.com/SnuJun/Kuggle)

#### Graduate
- [High-dimensional Statistical Inference](https://github.com/SnuJun/HighDim-Study)


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
