document.addEventListener("DOMContentLoaded", function () {
  // Create a lightbox overlay (single shared factory)
  function createLightbox() {
    if (document.querySelector('.wahi-lightbox')) return null;

    const overlay = document.createElement('div');
    overlay.className = 'wahi-lightbox';
    Object.assign(overlay.style, {
      position: 'fixed',
      top: 0,
      left: 0,
      right: 0,
      bottom: 0,
      background: 'rgba(0,0,0,0.85)',
      display: 'flex',
      alignItems: 'center',
      justifyContent: 'center',
      zIndex: 9999,
      padding: '20px',
      boxSizing: 'border-box',
    });

    const imgWrap = document.createElement('div');
    Object.assign(imgWrap.style, {
      maxWidth: '90%',
      maxHeight: '90%',
      textAlign: 'center',
    });

    const img = document.createElement('img');
    Object.assign(img.style, {
      maxWidth: '100%',
      maxHeight: '80vh',
      boxShadow: '0 8px 30px rgba(0,0,0,0.5)',
      borderRadius: '4px',
    });

    const caption = document.createElement('div');
    Object.assign(caption.style, {
      color: '#fff',
      marginTop: '12px',
      fontSize: '16px',
      lineHeight: '1.3',
      textAlign: 'center',
    });

    imgWrap.appendChild(img);
    imgWrap.appendChild(caption);
    overlay.appendChild(imgWrap);

    function removeOverlay() {
      if (overlay.parentNode) overlay.parentNode.removeChild(overlay);
      document.removeEventListener('keydown', escListener);
    }

    overlay.addEventListener('click', (e) => {
      if (e.target === overlay) removeOverlay();
    });

    function escListener(e) {
      if (e.key === 'Escape') removeOverlay();
    }

    document.addEventListener('keydown', escListener);

    return { overlay, img, caption };
  }

  // Initialize all slider instances on the page
  const sliders = document.querySelectorAll('.wahi-slider');
  if (!sliders || sliders.length === 0) return;

  sliders.forEach((slider) => {
    const track = slider.querySelector('.wahi-slider-track');
    if (!track) return;

    const slides = track.querySelectorAll('.wahi-slide');
    if (!slides || slides.length === 0) return;

    const prev = slider.querySelector('.wahi-slider-btn.prev');
    const next = slider.querySelector('.wahi-slider-btn.next');

    let index = 0;

    function update() {
      track.style.transform = `translateX(-${index * 100}%)`;
    }

    if (prev) {
      prev.addEventListener('click', () => {
        index = (index - 1 + slides.length) % slides.length;
        update();
      });
    }

    if (next) {
      next.addEventListener('click', () => {
        index = (index + 1) % slides.length;
        update();
      });
    }

    // attach click-to-enlarge for this slider's slides
    slides.forEach((slide) => {
      const image = slide.querySelector('img');
      const figcap = slide.querySelector('figcaption');
      if (!image) return;

      image.style.cursor = 'zoom-in';
      image.addEventListener('click', () => {
        const existing = document.querySelector('.wahi-lightbox');
        if (existing) return;

        const lb = createLightbox();
        if (!lb) return;
        lb.img.src = image.src;
        lb.img.alt = image.alt || '';
        lb.caption.textContent = figcap ? figcap.textContent : '';

        document.body.appendChild(lb.overlay);
      });
    });
  });
});