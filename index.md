<p align="center">
  <img src="/assets/img/profile.jpg" class="profile-pic" alt="Profile photo">
  <div class="profile-name">Junwoo Lee</div>
</p>
---
layout: default
title: Home
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
