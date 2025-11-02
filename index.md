# Welcome

Hello, I’m **Junwoo Lee**.  
I am currently pursuing a Master’s degree in Statistics at Seoul National University,  
where I study **high-dimensional statistics** and related topics in modern inference.

---

### Academic Background
- **Undergraduate**  
  - Regression Project  
  - *Economic Data Analysis*  
  - Machine Learning Club (*Kuggle*)
- **Graduate**  
  - High-dimensional Statistical Inference

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
