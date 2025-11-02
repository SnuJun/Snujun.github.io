<div style="display:flex; align-items:center; gap:20px;">
  <img src="/assets/img/profile.jpg" style="width:120px; height:120px; border-radius:50%; object-fit:cover; box-shadow:0 2px 6px rgba(0,0,0,0.2);">
  <div>
    <h1 style="margin-bottom:5px;">SnuJun’s Blog</h1>
    <p style="margin:0;">Statistics, ML, and life at SNU</p>
    <p style="margin-top:8px; font-weight:600;">Junwoo Lee</p>
  </div>
</div>

---

# 👋 Welcome!
Hi, I'm **Junwoo**. I write about statistics, high-dimensional inference, and things I’m learning.

- 👉 [About](/about)
- 👉 RSS: [/feed.xml](/feed.xml)

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
