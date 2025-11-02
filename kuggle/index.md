---
layout: default
title: Machine Learning Club Kuggle
permalink: /kuggle/
---

# Kuggle Sessions (Weeks 1–9)

<p>Click a week to open the PDF.</p>

<ul>
{% assign files = site.static_files | where_exp: "f", "f.path contains '/docs/'" %}
{% assign kuggle = files | where_exp: "f", "f.name contains 'Kuggle_wk' and f.extname == '.pdf'" %}
{% assign sorted = kuggle | sort: "name" %}
{% for f in sorted %}
  <li><a href="{{ f.path | relative_url }}" target="_blank" rel="noopener">
    {{ f.name | replace: 'Kuggle_', '' | replace: '_ppt', '' | replace: '.pdf','' | replace: '_', ' ' | capitalize }}
  </a></li>
{% endfor %}
</ul>
