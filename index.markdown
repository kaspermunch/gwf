---
layout: default
---

GWF (Grid WorkFlow) is a small utility for specifying and running workflows on a computer grid. Workflows are specified in Python, with some special functions to specify tasks to run, and GWF then keeps track of taskts that needs to be done in a make-like fashion and schedules jobs when they need to be run.


<div class="home">

  <h1 class="page-heading">Posts</h1>

  <ul class="post-list">
    {% for post in site.posts %}
      <li>
        <span class="post-meta">{{ post.date | date: "%b %-d, %Y" }}</span>

        <h2>
          <a class="post-link" href="{{ post.url | prepend: site.baseurl }}">{{ post.title }}</a>
        </h2>
      </li>
    {% endfor %}
  </ul>

  <p class="rss-subscribe">subscribe <a href="{{ "/feed.xml" | prepend: site.baseurl }}">via RSS</a></p>

</div>
