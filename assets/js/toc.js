// this file was created as an ad-hoc solution to a problem stemming from an inability to rebuild the scripts.min.js file via Grunt. It circumvents the
// minified js file and connects directly with the _toc.html file in the _includes directory. while this
// probably isn't a 'worst practice', it's definitely not a best practice, so tread lightly down similar paths!


$(document).ready(function() {
    var screenHeight = window.innerHeight || document.documentElement.clientHeight || document.body.clientHeight;

    window.addEventListener('scroll', function(event) {
      if (!$('#markdown-toc').length)
        return;

      var min_offset = null;
      var target_id = null;

      $(':header').each(function (index, item) {
        if ($(item).attr('id')) {
          var offset = $(item).offset().top - $(document).scrollTop();

          if (offset < (screenHeight / 2)) {
            if (min_offset == null || Math.abs(offset) < Math.abs(min_offset)) {
              min_offset = offset;
              target_id = $(item).attr('id');
            }
          }
        };
      });

      if (target_id && target_id != 'toc-header') {
        $("#markdown-toc a").removeClass('header-active');
        $("#markdown-toc a[href$='#"+target_id+"']").addClass('header-active');
      }
    }, true);

    // this is for smooth scrolling
    $('#markdown-toc a[href^="#"]').on('click', function(event) {
        var anchor = this;
        var target = $(anchor.getAttribute('href'));
        if( target.length ) {
            event.preventDefault();
            $('html, body').stop().animate({
                scrollTop: target.offset().top
            }, 500)
            .promise()
            .done(function(){ location.href = anchor.getAttribute('href'); });
        }
    });
  });