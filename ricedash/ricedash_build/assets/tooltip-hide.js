// assets/zz_tooltip-hide.js
(function () {
    const HIDE_DELAY_MS = 140;
    const IDLE_HIDE_MS  = 1200;
  
    function bind(container, tip) {
      if (!container || !tip) return;
      if (container.__tooltipHideBound) return;
      container.__tooltipHideBound = true;
  
      let hideTimer = null;
      let lastPointerTick = Date.now();
  
      function getNodeCanvas() {
        const byData = container.querySelector('canvas[data-id*="node"]');
        if (byData) return byData;
        const all = container.querySelectorAll('canvas');
        return all.length ? all[all.length - 1] : null;
      }
  
      let nodeCanvas = getNodeCanvas();
      let nodeCtx = nodeCanvas ? nodeCanvas.getContext('2d', { willReadFrequently: true }) : null;
  
      function refreshNodeCanvasIfNeeded() {
        const nc = getNodeCanvas();
        if (nc !== nodeCanvas) {
          nodeCanvas = nc;
          nodeCtx = nodeCanvas ? nodeCanvas.getContext('2d', { willReadFrequently: true }) : null;
        }
      }
  
      function isOverNode(e) {
        refreshNodeCanvasIfNeeded();
        if (!nodeCanvas || !nodeCtx) return false;
  
        const rect = nodeCanvas.getBoundingClientRect();
        const scaleX = nodeCanvas.width / rect.width;
        const scaleY = nodeCanvas.height / rect.height;
        const x = Math.floor((e.clientX - rect.left) * scaleX);
        const y = Math.floor((e.clientY - rect.top) * scaleY);
        if (x < 0 || y < 0 || x >= nodeCanvas.width || y >= nodeCanvas.height) return false;
  
        const a = nodeCtx.getImageData(x, y, 1, 1).data[3];
        return a > 0;
      }
  
      function cancelHide() {
        if (hideTimer) { clearTimeout(hideTimer); hideTimer = null; }
      }
  
      function scheduleHide() {
        cancelHide();
        hideTimer = setTimeout(() => {
          tip.style.display = 'none';
          tip.style.visibility = 'hidden';
          tip.style.opacity = '0';
        }, HIDE_DELAY_MS);
      }
  
      function onMoveLike(e) {
        lastPointerTick = Date.now();
        if (isOverNode(e)) {
          cancelHide(); // 노드 위면 유지
        } else {
          scheduleHide(); // 배경이면 숨김 예약
        }
      }
  
      container.addEventListener('mousemove', onMoveLike);
      container.addEventListener('touchmove', (ev) => {
        const t = ev.touches && ev.touches[0];
        if (!t) return;
        onMoveLike({ clientX: t.clientX, clientY: t.clientY });
      }, { passive: true });
  
      container.addEventListener('mouseleave', () => {
        tip.style.display = 'none';
        tip.style.visibility = 'hidden';
        tip.style.opacity = '0';
      });
  
      setInterval(() => {
        if (Date.now() - lastPointerTick > IDLE_HIDE_MS) {
          tip.style.display = 'none';
          tip.style.visibility = 'hidden';
          tip.style.opacity = '0';
        }
      }, 300);
  
      const obs = new MutationObserver(() => refreshNodeCanvasIfNeeded());
      obs.observe(container, { childList: true, subtree: true });
    }
  
    function tryBind() {
      const container = document.getElementById('pedigree-cytoscape');
      const tip = document.getElementById('tooltip');
      if (container && tip) bind(container, tip);
    }
  
    if (document.readyState === 'loading') {
      document.addEventListener('DOMContentLoaded', tryBind);
    } else {
      tryBind();
    }
  
    const rootObs = new MutationObserver(() => tryBind());
    rootObs.observe(document.body, { childList: true, subtree: true });
  })();
  